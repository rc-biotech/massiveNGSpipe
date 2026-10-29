/*
 * massiveNGSpipe trimcollapse: fastp-compatible adapter/quality trimming fused with read collapsing.
 *
 * One pass over a (gzipped) FASTQ: per read, in fastp's order
 *   1. trim front/tail          (fastp --trim_front1 / --trim_tail1)
 *   2. polyG tail trimming       (fastp default for two-colour instruments: auto | on | off)
 *   3. adapter trimming          (fastp AdapterTrimmer::trimBySequence semantics)
 *   4. optional tail trim after adapter removal (massiveNGSpipe 3' barcode/UMI case)
 *   5. filtering                 (fastp defaults: qualified phred 15, <= 40% unqualified, <= 5 N, length_required)
 * passing reads are counted in a hash and written as an ORFik collapsed fasta (>seq<i>_x<count>, by count), and a
 * fastp-style JSON report is written (summary, filtering_result, duplication, adapter_cutting,
 * read1_before/after_filtering with quality/content curves, command). An optional FASTQ of trimmed reads
 * (e.g. the first N reads, for barcode detection) can be written instead of / besides the collapsed fasta.
 *
 * Threads: one reader (decompression with ISA-L igzip when built with -DHAVE_ISAL, else zlib) cuts record-aligned blocks;
 * -w workers parse, compute statistics, trim and filter; the main thread collapses. Results do not depend on -w.
 *
 * Build: cc -O3 -o trimcollapse trimcollapse.c -lz -lpthread [-DHAVE_ISAL -Iisal -l:libisal.so.2]
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <time.h>
#include <pthread.h>
#ifdef HAVE_ISAL
#include "igzip_lib.h"
#endif

#define VERSION "0.1"
#define MAXLEN 4096
static void die(const char *m) { fprintf(stderr, "trimcollapse: %s\n", m); exit(2); }
static void *xmalloc(size_t n) { void *p = malloc(n); if (!p) die("out of memory"); return p; }

/* ---------------- options ---------------- */
typedef struct {
  const char *in, *collapsed, *json, *fastq_out, *adapter;
  int front, tail, tail_after_adapter, dimer_max, qf, minlen, qual_phred, unqual_pct, n_limit, polyg, polyg_minlen;
  long long max_reads;
  int workers;
  char cmd[8192];
} opts_t;

/* ---------------- statistics (fastp-style) ---------------- */
typedef struct {
  uint64_t reads, bases, q20, q30, q40, gc, maxlen;
  uint64_t *cnt[5];      /* per position base counts A,T,C,G,N */
  uint64_t *qsum[5];     /* per position quality sums by base */
  uint64_t *qall;        /* per position quality sum, all bases */
  uint64_t *nall;        /* per position base count, all */
  uint64_t kmer[1024];   /* 5-mers */
} stats_t;
static void stats_init(stats_t *s) {
  memset(s, 0, sizeof *s);
  for (int b = 0; b < 5; b++) { s->cnt[b] = calloc(MAXLEN, 8); s->qsum[b] = calloc(MAXLEN, 8); }
  s->qall = calloc(MAXLEN, 8); s->nall = calloc(MAXLEN, 8);
}
static int bidx(char c) { switch (c) { case 'A': return 0; case 'T': return 1; case 'C': return 2; case 'G': return 3; default: return 4; } }
static void stats_add(stats_t *s, const char *seq, const char *q, int len) {
  s->reads++; s->bases += len; if ((uint64_t)len > s->maxlen) s->maxlen = len;
  unsigned k = 0; int kl = 0;
  for (int i = 0; i < len; i++) {
    int b = bidx(seq[i]); int qq = q[i] - 33;
    s->cnt[b][i]++; s->qsum[b][i] += qq; s->qall[i] += qq; s->nall[i]++;
    if (qq >= 20) s->q20++; if (qq >= 30) s->q30++; if (qq >= 40) s->q40++;
    if (b == 2 || b == 3) s->gc++;
    if (b < 4) { k = ((k << 2) | (unsigned)b) & 1023; if (++kl >= 5) s->kmer[k]++; } else kl = 0;
  }
}

/* ---------------- collapse hash ---------------- */
typedef struct { char *s; uint32_t len; uint64_t n; uint64_t h; } entry;
typedef struct { entry *T; size_t cap, n; } table;            /* one per collapse shard, owned by one thread */
static uint64_t hash_bytes(const char *s, size_t n) { uint64_t h = 1469598103934665603ULL; for (size_t i = 0; i < n; i++) { h ^= (unsigned char)s[i]; h *= 1099511628211ULL; } return h ^ (h >> 29); }
static void tgrow(table *t) {
  size_t oc = t->cap; entry *O = t->T; t->cap = oc ? oc * 2 : (1 << 20); t->T = calloc(t->cap, sizeof(entry)); if (!t->T) die("out of memory");
  for (size_t i = 0; i < oc; i++) if (O[i].s) { size_t j = (O[i].h >> 8) & (t->cap - 1); while (t->T[j].s) j = (j + 1) & (t->cap - 1); t->T[j] = O[i]; }
  free(O);
}
static void tadd(table *t, const char *s, uint32_t len, uint64_t h) {
  if (t->n * 2 >= t->cap) tgrow(t);
  size_t i = (h >> 8) & (t->cap - 1);
  while (t->T[i].s) { if (t->T[i].h == h && t->T[i].len == len && !memcmp(t->T[i].s, s, len)) { t->T[i].n++; return; } i = (i + 1) & (t->cap - 1); }
  t->T[i].s = xmalloc(len); memcpy(t->T[i].s, s, len); t->T[i].len = len; t->T[i].n = 1; t->T[i].h = h; t->n++;
}
static int cmp_count(const void *a, const void *b) { const entry *x = a, *y = b; return (x->n < y->n) - (x->n > y->n); }

/* ---------------- fastp algorithms ---------------- */
/* PolyX::trimPolyG (fastp): scan from the 3' end; one mismatch allowed per 8 bases, at most 5 */
static int trim_polyg(const char *s, int len, int minlen) {
  const int one_mm_per = 8, max_mm = 5; int mm = 0, i, first_g = len - 1;
  for (i = 0; i < len; i++) {
    if (s[len - i - 1] != 'G') mm++; else first_g = len - i - 1;
    int allowed = (i + 1) / one_mm_per;
    if (mm > max_mm || (mm > allowed && i >= minlen - 1)) break;
  }
  if (i >= minlen) return first_g;     /* new length */
  return len;
}
/* Matcher::matchWithOneInsertion (fastp v1.3.7, verbatim logic): insData has one extra base somewhere within cmplen */
static int match_one_insertion(const char *ins, const char *norm, int cmplen, int diff_limit) {
  if (cmplen <= 0) return 0;
  int left[cmplen], right[cmplen];
  left[0] = ins[0] == norm[0] ? 0 : 1;
  right[cmplen - 1] = ins[cmplen] == norm[cmplen - 1] ? 0 : 1;
  for (int i = 1; i < cmplen; i++) {
    left[i] = left[i - 1] + (ins[i] != norm[i]);
    if (left[i] + right[cmplen - 1] > diff_limit) break;
  }
  for (int i = cmplen - 2; i >= 0; i--) {
    right[i] = right[i + 1] + (ins[i + 1] != norm[i]);
    if (right[i] + left[0] > diff_limit) { for (int q = 0; q < i; q++) right[q] = diff_limit + 1; break; }
  }
  for (int i = 1; i < cmplen; i++) {
    if (left[i - 1] + right[cmplen - 1] > diff_limit) return 0;
    if (left[i - 1] + right[i] <= diff_limit) return 1;
  }
  return 0;
}
/* AdapterTrimmer::trimBySequence (fastp v1.3.7): (1) hamming match from a small negative offset (skipped A-tail), one
   mismatch allowed per 8 compared bases; (2) if none, one insertion, then (3) one deletion -- as in fastp, steps 2/3
   compare the START of the read with the adapter for every pos (fastp passes rdata, not rdata + pos), and cut at that
   pos. Reproduced as is so reports and outputs match fastp 1.3.7. Returns the new length (0 = whole read), -1 if none. */
static int trim_adapter(const char *s, int len, const char *ad, int alen) {
  const int match_req = 4, one_mm_per = 8;
  if (alen < match_req) return -1;
  int start = alen >= 16 ? -4 : alen >= 12 ? -3 : alen >= 8 ? -2 : 0;
  for (int pos = start; pos < len - match_req; pos++) {
    int cmplen = len - pos < alen ? len - pos : alen;
    int allowed = cmplen / one_mm_per, mm = 0, ok = 1;
    for (int i = pos < 0 ? -pos : 0; i < cmplen; i++) if (ad[i] != s[i + pos]) { if (++mm > allowed) { ok = 0; break; } }
    if (ok) return pos < 0 ? 0 : pos;
  }
  for (int pos = 0; pos < len - match_req - 1; pos++) {           /* one insertion in the read */
    int cmplen = len - pos - 1 < alen ? len - pos - 1 : alen;
    if (match_one_insertion(s, ad, cmplen, cmplen / one_mm_per - 1)) return pos;
  }
  for (int pos = 0; pos < len - match_req; pos++) {               /* one deletion in the read */
    int cmplen = len - pos < alen - 1 ? len - pos : alen - 1;
    if (match_one_insertion(ad, s, cmplen, cmplen / one_mm_per - 1)) return pos;
  }
  return -1;
}

/* ---------------- adapter sequences actually cut (fastp read1_adapter_counts) ---------------- */
typedef struct { char *s; uint64_t n; } acount;
static acount *AC; static size_t acn, accap;
static void adapter_seen(const char *s, int len) {
  if (len > 48) len = 48; if (len < 1) return;
  for (size_t i = 0; i < acn; i++) if ((int)strlen(AC[i].s) == len && !memcmp(AC[i].s, s, len)) { AC[i].n++; return; }
  if (acn < 100000) { if (acn == accap) { accap = accap ? accap * 2 : 1024; AC = realloc(AC, accap * sizeof(acount)); } AC[acn].s = strndup(s, len); AC[acn].n = 1; acn++; }
}

/* ---------------- JSON ---------------- */
static void json_curves(FILE *f, stats_t *s) {
  const char *names[5] = {"A", "T", "C", "G", "N"}; int L = s->maxlen;
  fprintf(f, "\t\t\"quality_curves\": {\n");
  for (int b = 0; b < 4; b++) {
    fprintf(f, "\t\t\t\"%s\":[", names[b]);
    for (int i = 0; i < L; i++) fprintf(f, "%s%g", i ? "," : "", s->cnt[b][i] ? (double)s->qsum[b][i] / s->cnt[b][i] : 0.0);
    fprintf(f, "],\n");
  }
  fprintf(f, "\t\t\t\"mean\":[");
  for (int i = 0; i < L; i++) fprintf(f, "%s%g", i ? "," : "", s->nall[i] ? (double)s->qall[i] / s->nall[i] : 0.0);
  fprintf(f, "]\n\t\t},\n\t\t\"content_curves\": {\n");
  for (int b = 0; b < 5; b++) {
    fprintf(f, "\t\t\t\"%s\":[", names[b]);
    for (int i = 0; i < L; i++) fprintf(f, "%s%g", i ? "," : "", s->nall[i] ? (double)s->cnt[b][i] / s->nall[i] : 0.0);
    fprintf(f, "],\n");
  }
  fprintf(f, "\t\t\t\"GC\":[");
  for (int i = 0; i < L; i++) fprintf(f, "%s%g", i ? "," : "", s->nall[i] ? (double)(s->cnt[2][i] + s->cnt[3][i]) / s->nall[i] : 0.0);
  fprintf(f, "]\n\t\t},\n");
}
static void json_read_block(FILE *f, const char *name, stats_t *s, int last) {
  const char nt[] = "ATCG";
  fprintf(f, "\t\"%s\": {\n\t\t\"total_reads\":%llu,\n\t\t\"total_bases\":%llu,\n\t\t\"q20_bases\":%llu,\n\t\t\"q30_bases\":%llu,\n\t\t\"q40_bases\":%llu,\n\t\t\"total_cycles\":%llu,\n",
          name, (unsigned long long)s->reads, (unsigned long long)s->bases, (unsigned long long)s->q20, (unsigned long long)s->q30,
          (unsigned long long)s->q40, (unsigned long long)s->maxlen);
  json_curves(f, s);
  fprintf(f, "\t\t\"kmer_count\": {\n");
  for (int k = 0; k < 1024; k++) {
    char km[6]; for (int j = 0; j < 5; j++) km[j] = nt[(k >> (2 * (4 - j))) & 3]; km[5] = 0;
    fprintf(f, "\t\t\t\"%s\":%llu%s\n", km, (unsigned long long)s->kmer[k], k < 1023 ? "," : "");
  }
  fprintf(f, "\t\t},\n\t\t\"overrepresented_sequences\": {\n\t\t}\n\t}%s\n", last ? "" : ",");
}
static void json_summary(FILE *f, const char *name, stats_t *s, int last) {
  fprintf(f, "\t\t\"%s\": {\n\t\t\t\"total_reads\":%llu,\n\t\t\t\"total_bases\":%llu,\n\t\t\t\"q20_bases\":%llu,\n\t\t\t\"q30_bases\":%llu,\n"
             "\t\t\t\"q20_rate\":%g,\n\t\t\t\"q30_rate\":%g,\n\t\t\t\"read1_mean_length\":%llu,\n\t\t\t\"gc_content\":%g\n\t\t}%s\n",
          name, (unsigned long long)s->reads, (unsigned long long)s->bases, (unsigned long long)s->q20, (unsigned long long)s->q30,
          s->bases ? (double)s->q20 / s->bases : 0, s->bases ? (double)s->q30 / s->bases : 0,
          (unsigned long long)(s->reads ? s->bases / s->reads : 0), s->bases ? (double)s->gc / s->bases : 0, last ? "" : ",");
}
static int cmp_ac(const void *a, const void *b) { const acount *x = a, *y = b; return (x->n < y->n) - (x->n > y->n); }

/* ---------------- block queue ---------------- */
typedef struct block { char *buf; size_t len; struct block *next; } block;
typedef struct { block *head, *tail; int n, max, closed; pthread_mutex_t m; pthread_cond_t c_put, c_get; } queue;
static void q_init(queue *q, int max) { memset(q, 0, sizeof *q); q->max = max; pthread_mutex_init(&q->m, NULL); pthread_cond_init(&q->c_put, NULL); pthread_cond_init(&q->c_get, NULL); }
static void q_put(queue *q, block *b) {
  pthread_mutex_lock(&q->m); while (q->n >= q->max) pthread_cond_wait(&q->c_put, &q->m);
  b->next = NULL; if (q->tail) q->tail->next = b; else q->head = b; q->tail = b; q->n++;
  pthread_cond_signal(&q->c_get); pthread_mutex_unlock(&q->m);
}
static block *q_get(queue *q) {
  pthread_mutex_lock(&q->m); while (!q->n && !q->closed) pthread_cond_wait(&q->c_get, &q->m);
  block *b = q->head; if (b) { q->head = b->next; if (!q->head) q->tail = NULL; q->n--; pthread_cond_signal(&q->c_put); }
  pthread_mutex_unlock(&q->m); return b;
}
static void q_close(queue *q) { pthread_mutex_lock(&q->m); q->closed = 1; pthread_cond_broadcast(&q->c_get); pthread_mutex_unlock(&q->m); }

/* ---------------- reader: decompress, cut record-aligned blocks ---------------- */
#define BLOCK (8 << 20)
static opts_t O; static queue QIN, QOUT; static int NSHARD; static queue *QS; static table *TS; static int polyg_on; static char first_header[MAXLEN]; static volatile int reader_err;
typedef struct { FILE *f; int gz; unsigned char *ib; size_t ilen, ipos; int eof;
#ifdef HAVE_ISAL
  struct inflate_state st;
#else
  z_stream zs;
#endif
} src_t;
static size_t src_read(src_t *r, char *out, size_t cap) {      /* returns bytes of decompressed data, 0 at end */
  if (!r->gz) return fread(out, 1, cap, r->f);
  size_t got = 0;
  while (got < cap) {
    if (r->ipos >= r->ilen && !r->eof) { r->ilen = fread(r->ib, 1, 1 << 22, r->f); r->ipos = 0; if (!r->ilen) r->eof = 1; }
#ifdef HAVE_ISAL
    if (r->st.block_state == ISAL_BLOCK_FINISH) {                /* next gzip member, if any */
      if (r->ipos >= r->ilen && r->eof) break;
      isal_inflate_reset(&r->st); r->st.crc_flag = ISAL_GZIP;
    }
    r->st.next_in = r->ib + r->ipos; r->st.avail_in = r->ilen - r->ipos;
    r->st.next_out = (uint8_t *)out + got; r->st.avail_out = cap - got;
    int ret = isal_inflate(&r->st); if (ret < 0) die("gzip decompression failed (isa-l)");
    size_t used = (r->ilen - r->ipos) - r->st.avail_in; r->ipos += used;
    size_t made = (cap - got) - r->st.avail_out; got += made;
    if (!made && !used && r->eof && r->st.block_state != ISAL_BLOCK_FINISH) break;
#else
    r->zs.next_in = r->ib + r->ipos; r->zs.avail_in = r->ilen - r->ipos;
    r->zs.next_out = (Bytef *)out + got; r->zs.avail_out = cap - got;
    int ret = inflate(&r->zs, Z_NO_FLUSH);
    size_t used = (r->ilen - r->ipos) - r->zs.avail_in; r->ipos += used;
    got += (cap - got) - r->zs.avail_out;
    if (ret == Z_STREAM_END) { if (r->ipos >= r->ilen && r->eof) break; inflateReset(&r->zs); }
    else if (ret != Z_OK && ret != Z_BUF_ERROR) die("gzip decompression failed (zlib)");
    if (!used && r->eof && r->zs.avail_in == 0 && ret == Z_BUF_ERROR) break;
#endif
  }
  return got;
}
static void *reader(void *arg) {
  (void)arg; src_t r = {0}; r.f = fopen(O.in, "rb"); if (!r.f) die("cannot open input");
  size_t n = strlen(O.in); r.gz = n > 3 && !strcmp(O.in + n - 3, ".gz");
  if (r.gz) { r.ib = xmalloc(1 << 22);
#ifdef HAVE_ISAL
    isal_inflate_init(&r.st); r.st.crc_flag = ISAL_GZIP;
#else
    if (inflateInit2(&r.zs, 15 + 32) != Z_OK) die("zlib init");
#endif
  }
  char *carry = NULL; size_t clen = 0; long long reads = 0; int lines_mod = 0, done = 0;
  while (!done) {
    char *buf = xmalloc(BLOCK + MAXLEN * 4); size_t len = 0;
    if (clen) { memcpy(buf, carry, clen); len = clen; free(carry); carry = NULL; clen = 0; }
    size_t got = src_read(&r, buf + len, BLOCK - len); len += got;
    if (!got) done = 1;
    /* cut after the last complete record; count records (for -r) */
    size_t cut = 0; int lm = 0; long long rec = 0;
    for (char *q = buf, *e = buf + len; q < e; ) {
      char *nl = memchr(q, '\n', e - q); if (!nl) break; q = nl + 1;
      if (++lm == 4) { lm = 0; rec++; cut = q - buf; if (O.max_reads >= 0 && reads + rec >= O.max_reads) { done = 1; break; } }
    }
    (void)lines_mod;
    if (cut < len && !done) { clen = len - cut; carry = xmalloc(clen); memcpy(carry, buf + cut, clen); }
    else if (cut < len && done && O.max_reads < 0 && len - cut > 1) die("truncated FASTQ (incomplete last record)");
    reads += rec;
    if (!first_header[0] && cut) { size_t e = 0; while (e < cut && buf[e] != '\n' && e < MAXLEN - 1) e++; memcpy(first_header, buf, e); first_header[e] = 0; }
    if (!cut) { free(buf); continue; }
    block *b = xmalloc(sizeof(block)); b->buf = buf; b->len = cut; q_put(&QIN, b);
  }
  free(carry); fclose(r.f); q_close(&QIN); return NULL;
}

/* ---------------- workers ---------------- */
typedef struct { stats_t before, after; uint64_t lowq, nn, shortr, dimer, ad_reads, ad_bases; acount *ac; size_t acn, accap; } wstate;
static void w_adapter_seen(wstate *w, const char *s, int len) {
  if (len > 48) len = 48; if (len < 1) return;
  for (size_t i = 0; i < w->acn && i < 64; i++) if ((int)strlen(w->ac[i].s) == len && !memcmp(w->ac[i].s, s, len)) { w->ac[i].n++; return; }
  if (w->acn < 4096) { if (w->acn == w->accap) { w->accap = w->accap ? w->accap * 2 : 256; w->ac = realloc(w->ac, w->accap * sizeof(acount)); } w->ac[w->acn].s = strndup(s, len); w->ac[w->acn].n = 1; w->acn++; }
}
static void *worker(void *arg) {
  wstate *w = arg; int use_adapter = strcmp(O.adapter, "disable") != 0; int alen = use_adapter ? strlen(O.adapter) : 0;
  char qmin = (char)(O.qual_phred + 33);
  block *b;
  while ((b = q_get(&QIN))) {
    /* output: [uint16 len][seq] for collapse, or fastq records when --fastq-out */
    char *out = NULL; size_t olen = 0; char *p = b->buf, *end = b->buf + b->len;
    char **sb = NULL; size_t *sl = NULL, *sc = NULL;
    if (O.fastq_out) out = xmalloc(b->len + 16);
    else { sb = xmalloc(NSHARD * sizeof(char *)); sl = calloc(NSHARD, sizeof(size_t)); sc = xmalloc(NSHARD * sizeof(size_t));
           for (int k = 0; k < NSHARD; k++) { sc[k] = b->len / NSHARD * 2 + 65536; sb[k] = xmalloc(sc[k]); } }
    while (p < end) {
      char *h = p; char *s = memchr(h, '\n', end - h) + 1; char *pl = memchr(s, '\n', end - s) + 1;
      char *q = memchr(pl, '\n', end - pl) + 1; char *nx = memchr(q, '\n', end - q) + 1;
      int len = pl - s - 1; if (len && s[len - 1] == '\r') len--;
      int ql = nx - q - 1; if (ql && q[ql - 1] == '\r') ql--;
      if (ql != len) die("sequence and quality lengths differ");
      if (len >= MAXLEN) die("read longer than 4095 nt");
      p = nx;
      stats_add(&w->before, s, q, len);
      int st = O.front, en = len - O.tail; if (en < st) en = st;
      const char *rs = s + st, *rq = q + st; int rl = en - st;
      if (polyg_on && rl > 0) rl = trim_polyg(rs, rl, O.polyg_minlen);
      int cut = 0;
      if (use_adapter && rl > 0) { int nl = trim_adapter(rs, rl, O.adapter, alen); if (nl >= 0) { w->ad_reads++; w->ad_bases += rl - nl; w_adapter_seen(w, rs + nl, rl - nl); rl = nl; cut = 1; } }
      if (cut && rl <= O.dimer_max) { w->dimer++; continue; }
      if (rl <= 0) { w->shortr++; continue; }
      if (O.qf) {                                   /* fastp quality filter (off with -Q, as in the pipeline's main pass) */
        int lowq = 0, nn = 0;
        for (int i = 0; i < rl; i++) { lowq += rq[i] < qmin; nn += rs[i] == 'N'; }
        if (lowq > O.unqual_pct * rl / 100.0) { w->lowq++; continue; }
        if (nn > O.n_limit) { w->nn++; continue; }
      }
      if (rl < O.minlen) { w->shortr++; continue; }
      if (O.tail_after_adapter > 0) {               /* second fastp stage of fastp_post_adapter_tail(): tail cut, length. That stage
                                                       reads stdin, so fastp cannot detect the instrument and polyG trimming is off */
        rl = rl > O.tail_after_adapter ? rl - O.tail_after_adapter : 0;
        if (rl <= 0 || rl < O.minlen) { w->shortr++; continue; }
      }
      stats_add(&w->after, rs, rq, rl);
      if (O.fastq_out) { int hl = s - h; memcpy(out + olen, h, hl); olen += hl; memcpy(out + olen, rs, rl); olen += rl; memcpy(out + olen, "\n+\n", 3); olen += 3; memcpy(out + olen, rq, rl); olen += rl; out[olen++] = '\n'; }
      else {
        uint64_t hh = hash_bytes(rs, rl); int k = hh % NSHARD; uint16_t l16 = rl;
        /* grow: all copies of one sequence go to one shard, so a dominant sequence overflows a fixed share of the block */
        if (sl[k] + 10 + rl > sc[k]) { sc[k] = sc[k] * 2 + 10 + rl; sb[k] = realloc(sb[k], sc[k]); if (!sb[k]) die("out of memory"); }
        char *d = sb[k] + sl[k];
        memcpy(d, &hh, 8); memcpy(d + 8, &l16, 2); memcpy(d + 10, rs, rl); sl[k] += 10 + rl;
      }
    }
    free(b->buf);
    if (O.fastq_out) { b->buf = out; b->len = olen; q_put(&QOUT, b); }
    else { for (int k = 0; k < NSHARD; k++) { block *x = xmalloc(sizeof(block)); x->buf = sb[k]; x->len = sl[k]; q_put(&QS[k], x); } free(sb); free(sl); free(sc); free(b); }
  }
  return NULL;
}
static void stats_merge(stats_t *a, stats_t *b) {
  a->reads += b->reads; a->bases += b->bases; a->q20 += b->q20; a->q30 += b->q30; a->q40 += b->q40; a->gc += b->gc;
  if (b->maxlen > a->maxlen) a->maxlen = b->maxlen;
  for (int i = 0; i < MAXLEN; i++) { for (int k = 0; k < 5; k++) { a->cnt[k][i] += b->cnt[k][i]; a->qsum[k][i] += b->qsum[k][i]; } a->qall[i] += b->qall[i]; a->nall[i] += b->nall[i]; }
  for (int k = 0; k < 1024; k++) a->kmer[k] += b->kmer[k];
}

/* joins the workers, then closes the output queue so the collapsing loop ends */
static pthread_t *CLOSE_WT; static int CLOSE_N;
static void *close_fn(void *x) { (void)x; for (int i = 0; i < CLOSE_N; i++) pthread_join(CLOSE_WT[i], NULL); q_close(&QOUT);
  for (int k = 0; k < NSHARD; k++) q_close(&QS[k]); return NULL; }
static void *collector(void *x) {            /* one per shard: owns TS[k] */
  int k = (int)(intptr_t)x; block *b;
  while ((b = q_get(&QS[k]))) {
    char *p = b->buf, *e = b->buf + b->len;
    while (p < e) { uint64_t h; uint16_t l; memcpy(&h, p, 8); memcpy(&l, p + 8, 2); tadd(&TS[k], p + 10, l, h); p += 10 + l; }
    free(b->buf); free(b);
  }
  return NULL;
}

/* ---------------- main ---------------- */
int main(int argc, char **argv) {
  opts_t o = {0}; o.adapter = "disable"; o.minlen = 15; o.qual_phred = 15; o.unqual_pct = 40; o.n_limit = 5; o.polyg = -1; o.polyg_minlen = 10;
  o.max_reads = -1; o.dimer_max = 2; o.workers = 4; o.qf = 1;
  size_t cl = 0; for (int i = 0; i < argc; i++) cl += snprintf(o.cmd + cl, sizeof o.cmd - cl, "%s%s", i ? " " : "", argv[i]);
  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!strcmp(a, "-i")) o.in = argv[++i]; else if (!strcmp(a, "-o")) o.collapsed = argv[++i];
    else if (!strcmp(a, "-j")) o.json = argv[++i]; else if (!strcmp(a, "--fastq-out")) o.fastq_out = argv[++i];
    else if (!strcmp(a, "-a")) o.adapter = argv[++i]; else if (!strcmp(a, "-f")) o.front = atoi(argv[++i]);
    else if (!strcmp(a, "-t")) o.tail = atoi(argv[++i]); else if (!strcmp(a, "--tail-after-adapter")) o.tail_after_adapter = atoi(argv[++i]);
    else if (!strcmp(a, "-l")) o.minlen = atoi(argv[++i]); else if (!strcmp(a, "-q")) o.qual_phred = atoi(argv[++i]);
    else if (!strcmp(a, "-u")) o.unqual_pct = atoi(argv[++i]); else if (!strcmp(a, "-n")) o.n_limit = atoi(argv[++i]);
    else if (!strcmp(a, "--poly-g")) { const char *v = argv[++i]; o.polyg = !strcmp(v, "on") ? 1 : !strcmp(v, "off") ? 0 : -1; }
    else if (!strcmp(a, "-r")) o.max_reads = atoll(argv[++i]);
    else if (!strcmp(a, "--dimer-max-len")) o.dimer_max = atoi(argv[++i]);
    else if (!strcmp(a, "-w")) o.workers = atoi(argv[++i]);
    else if (!strcmp(a, "-Q")) o.qf = 0;
    else die("unknown argument");
  }
  if (!o.in || (!o.collapsed && !o.fastq_out)) die("usage: trimcollapse -i in.fastq[.gz] (-o collapsed.fasta.gz | --fastq-out trimmed.fastq) [-j report.json] [-a ADAPTER|disable] [-f N] [-t N] [--tail-after-adapter N] [-l MIN] [-r N] [-w workers]");
  if (o.workers < 1) o.workers = 1;
  O = o; int use_adapter = strcmp(o.adapter, "disable") != 0;
  /* polyG auto-detection needs the first header: peek it before the workers start */
  { char peek[MAXLEN] = {0}; gzFile g = gzopen(o.in, "rb"); if (!g || !gzgets(g, peek, MAXLEN)) die("cannot read input"); gzclose(g);
    polyg_on = o.polyg >= 0 ? o.polyg : (!strncmp(peek, "@NS", 3) || !strncmp(peek, "@NB", 3) || !strncmp(peek, "@NDX", 4) || !strncmp(peek, "@A0", 3) ||
                                         !strncmp(peek, "@VH", 3) || !strncmp(peek, "@LH", 3) || !strncmp(peek, "@VL", 3)); }
  q_init(&QIN, 2 * o.workers + 2); q_init(&QOUT, 2 * o.workers + 2);
  NSHARD = o.fastq_out ? 1 : (o.workers < 8 ? o.workers : 8);
  QS = calloc(NSHARD, sizeof(queue)); TS = calloc(NSHARD, sizeof(table)); pthread_t *ct = xmalloc(NSHARD * sizeof(pthread_t));
  for (int k = 0; k < NSHARD; k++) q_init(&QS[k], 4 * o.workers + 4);
  if (!o.fastq_out) for (int k = 0; k < NSHARD; k++) pthread_create(&ct[k], NULL, collector, (void *)(intptr_t)k);
  wstate *W = calloc(o.workers, sizeof(wstate));
  for (int i = 0; i < o.workers; i++) { stats_init(&W[i].before); stats_init(&W[i].after); }
  pthread_t rt, *wt = xmalloc(o.workers * sizeof(pthread_t));
  pthread_create(&rt, NULL, reader, NULL);
  for (int i = 0; i < o.workers; i++) pthread_create(&wt[i], NULL, worker, &W[i]);
  /* closer: when all workers finish, close the output queue */
  FILE *fq = o.fastq_out ? fopen(o.fastq_out, "w") : NULL; if (o.fastq_out && !fq) die("cannot write fastq");
  pthread_t closer; CLOSE_WT = wt; CLOSE_N = o.workers;
  pthread_create(&closer, NULL, close_fn, NULL);
  block *b;
  if (fq) while ((b = q_get(&QOUT))) { fwrite(b->buf, 1, b->len, fq); free(b->buf); free(b); }
  pthread_join(closer, NULL); pthread_join(rt, NULL); if (fq) fclose(fq);
  if (!o.fastq_out) for (int k = 0; k < NSHARD; k++) pthread_join(ct[k], NULL);
  size_t tn = 0; for (int k = 0; k < NSHARD; k++) tn += TS[k].n;
  stats_t before, after; stats_init(&before); stats_init(&after);
  uint64_t f_lowq = 0, f_n = 0, f_short = 0, f_dimer = 0, ad_reads = 0, ad_bases = 0;
  for (int i = 0; i < o.workers; i++) {
    stats_merge(&before, &W[i].before); stats_merge(&after, &W[i].after);
    f_lowq += W[i].lowq; f_n += W[i].nn; f_short += W[i].shortr; f_dimer += W[i].dimer; ad_reads += W[i].ad_reads; ad_bases += W[i].ad_bases;
    for (size_t k = 0; k < W[i].acn; k++) for (uint64_t r = 0; r < 1; r++) {
      int found = 0; for (size_t j = 0; j < acn; j++) if (!strcmp(AC[j].s, W[i].ac[k].s)) { AC[j].n += W[i].ac[k].n; found = 1; break; }
      if (!found) { if (acn == accap) { accap = accap ? accap * 2 : 1024; AC = realloc(AC, accap * sizeof(acount)); } AC[acn].s = W[i].ac[k].s; AC[acn].n = W[i].ac[k].n; acn++; }
    }
  }
  /* collapsed fasta */
  if (o.collapsed) {
    entry *E = xmalloc((tn ? tn : 1) * sizeof(entry)); size_t k = 0;
    for (int sh = 0; sh < NSHARD; sh++) for (size_t i = 0; i < TS[sh].cap; i++) if (TS[sh].T[i].s) E[k++] = TS[sh].T[i];
    qsort(E, k, sizeof(entry), cmp_count);
    char tmpf[4096]; snprintf(tmpf, sizeof tmpf, "%s.tmp", o.collapsed);
    char cmd[8300]; snprintf(cmd, sizeof cmd, "pigz -p %d -6 > '%s'", o.workers, tmpf);
    FILE *pz = system("command -v pigz >/dev/null 2>&1") == 0 ? popen(cmd, "w") : NULL;
    char *ob = xmalloc(1 << 22); size_t ol = 0;
    gzFile g = NULL; if (!pz) { g = gzopen(tmpf, "wb6"); if (!g) die("cannot write collapsed output"); gzbuffer(g, 1 << 22); }
    for (size_t i = 0; i < k; i++) {
      if (ol > (1 << 22) - 8192) { if (pz) fwrite(ob, 1, ol, pz); else gzwrite(g, ob, ol); ol = 0; }
      ol += sprintf(ob + ol, ">seq%zu_x%llu\n", i + 1, (unsigned long long)E[i].n); memcpy(ob + ol, E[i].s, E[i].len); ol += E[i].len; ob[ol++] = '\n';
    }
    if (pz) { fwrite(ob, 1, ol, pz); if (pclose(pz) != 0) die("pigz failed"); } else { gzwrite(g, ob, ol); if (gzclose(g) != Z_OK) die("gzclose failed"); }
    if (rename(tmpf, o.collapsed)) die("writing collapsed output failed");
  }
  /* json */
  if (o.json) {
    char tmpj[4096]; snprintf(tmpj, sizeof tmpj, "%s.tmp", o.json);
    FILE *f = fopen(tmpj, "w"); if (!f) die("cannot write json");
    fprintf(f, "{\n\t\"summary\": {\n\t\t\"fastp_version\": \"massiveNGSpipe-trimcollapse %s (fastp-compatible)\",\n\t\t\"sequencing\": \"single end (%llu cycles)\",\n",
            VERSION, (unsigned long long)before.maxlen);
    json_summary(f, "before_filtering", &before, 0); json_summary(f, "after_filtering", &after, 1);
    fprintf(f, "\t},\n\t\"filtering_result\": {\n\t\t\"passed_filter_reads\": %llu,\n\t\t\"low_quality_reads\": %llu,\n\t\t\"too_many_N_reads\": %llu,\n",
            (unsigned long long)after.reads, (unsigned long long)f_lowq, (unsigned long long)f_n);
    if (use_adapter) fprintf(f, "\t\t\"adapter_dimer_reads\": %llu,\n", (unsigned long long)f_dimer);   /* fastp omits it without adapter trimming */
    fprintf(f, "\t\t\"too_short_reads\": %llu,\n\t\t\"too_long_reads\": 0\n\t},\n", (unsigned long long)f_short);
    fprintf(f, "\t\"duplication\": {\n\t\t\"rate\": %g\n\t},\n", (o.collapsed && after.reads) ? 1.0 - (double)tn / after.reads : 0.0);
    if (use_adapter) {                    /* fastp writes adapter_cutting only when adapter trimming is on */
    fprintf(f, "\t\"adapter_cutting\": {\n\t\t\"adapter_trimmed_reads\": %llu,\n\t\t\"adapter_trimmed_bases\": %llu,\n", (unsigned long long)ad_reads, (unsigned long long)ad_bases);
    if (use_adapter) fprintf(f, "\t\t\"read1_adapter_sequence\": \"%s\",\n", o.adapter);
    qsort(AC, acn, sizeof(acount), cmp_ac); uint64_t other = 0;
    fprintf(f, "\t\t\"read1_adapter_counts\": {\n");
    for (size_t i = 0; i < acn; i++) { if (i < 10) fprintf(f, "\t\t\t\"%s\": %llu,\n", AC[i].s, (unsigned long long)AC[i].n); else other += AC[i].n; }
    fprintf(f, "\t\t\t\"others\": %llu\n\t\t}\n\t},\n", (unsigned long long)other);
    }
    json_read_block(f, "read1_before_filtering", &before, 0); json_read_block(f, "read1_after_filtering", &after, 0);
    fprintf(f, "\t\"command\": \"");
    for (char *c = o.cmd; *c; c++) { if (*c == '"' || *c == '\\') fputc('\\', f); fputc(*c, f); }
    fprintf(f, "\"\n}\n");
    fclose(f); if (rename(tmpj, o.json)) die("rename json failed");
  }
  fprintf(stderr, "trimcollapse: %llu reads, %llu passed (%.1f%%), adapter in %llu, unique %zu\n", (unsigned long long)before.reads,
          (unsigned long long)after.reads, before.reads ? 100.0 * after.reads / before.reads : 0, (unsigned long long)ad_reads, tn);
  return 0;
}
