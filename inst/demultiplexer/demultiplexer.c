/*
 * massiveNGSpipe demultiplexer: fast pooled Ribo-seq demultiplexing and read collapsing.
 *
 * Read layout (read 1):  [UMI5][insert][UMI3][BARCODE][ADAPTER...]
 *
 * demux    -i pool.fastq[.gz] -b barcodes.tsv -o outdir -a ADAPTER -u5 N -u3 N -m MIN [-t TAG]
 *          [-am adapter_mismatches] [-bm barcode_mismatches] [-ak adapter_kmer_length]
 *   barcodes.tsv: "<sample>\t<barcode>" per line (all barcodes the same length).
 *   For every read: find the first adapter hit (first -ak bases of ADAPTER, <= -am mismatches) at a position
 *   that leaves room for UMIs + barcode + MIN insert bases; take the barcode right before it (exact match, or a
 *   unique match within -bm mismatches); strip both UMIs; keep inserts of >= MIN bases.
 *   Writes <outdir>/<sample>.<TAG>.seq (one insert per line) and <outdir>/stats.<TAG>.tsv.
 *   Input is read through `pigz -dc` (or `gzip -dc`) when it ends in .gz.
 *
 * collapse -o out.fasta.gz -m MIN in1.seq [in2.seq ...]
 *   Counts identical inserts over all inputs and writes ORFik "ribotoolkit" collapsed fasta
 *   (>seq<i>_x<count>, sorted by decreasing count), gzip compressed.
 *
 * Build: cc -O3 -o demultiplexer demultiplexer.c -lz
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

static void die(const char *m) { fprintf(stderr, "demultiplexer: %s\n", m); exit(2); }
static void *xmalloc(size_t n) { void *p = malloc(n); if (!p) die("out of memory"); return p; }
static void *xrealloc(void *p, size_t n) { p = realloc(p, n); if (!p) die("out of memory"); return p; }

/* ---------- open-addressing hash: string -> int (barcode lookup) ---------- */
typedef struct { char **key; int *val; size_t cap; } smap;
static uint64_t hash_bytes(const char *s, size_t n) { uint64_t h = 1469598103934665603ULL; for (size_t i = 0; i < n; i++) { h ^= (unsigned char)s[i]; h *= 1099511628211ULL; } return h; }
static void smap_init(smap *m, size_t cap) { m->cap = cap; m->key = calloc(cap, sizeof(char *)); m->val = calloc(cap, sizeof(int)); if (!m->key || !m->val) die("out of memory"); }
static int smap_get(smap *m, const char *k, size_t n) {
  size_t i = hash_bytes(k, n) & (m->cap - 1);
  while (m->key[i]) { if (!memcmp(m->key[i], k, n)) return m->val[i]; i = (i + 1) & (m->cap - 1); }
  return -1;
}
static void smap_put(smap *m, const char *k, size_t n, int v) {
  size_t i = hash_bytes(k, n) & (m->cap - 1);
  while (m->key[i]) { if (!memcmp(m->key[i], k, n)) { m->val[i] = v; return; } i = (i + 1) & (m->cap - 1); }
  m->key[i] = xmalloc(n + 1); memcpy(m->key[i], k, n); m->key[i][n] = 0; m->val[i] = v;
}

/* ---------- demux ---------- */
typedef struct { char *name; FILE *out; uint64_t exact, rescued, shortn, len_sum, raw_len_sum; } sample_t;

static int run_demux(int argc, char **argv) {
  const char *in = NULL, *bcf = NULL, *outdir = NULL, *adapter = "AGATCGGAAGAGC", *tag = "pool";
  int u5 = 0, u3 = 0, minlen = 15, am = 2, bm = 1, ak = 12;
  for (int i = 2; i < argc; i++) {
    if (!strcmp(argv[i], "-i")) in = argv[++i]; else if (!strcmp(argv[i], "-b")) bcf = argv[++i];
    else if (!strcmp(argv[i], "-o")) outdir = argv[++i]; else if (!strcmp(argv[i], "-a")) adapter = argv[++i];
    else if (!strcmp(argv[i], "-t")) tag = argv[++i]; else if (!strcmp(argv[i], "-u5")) u5 = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-u3")) u3 = atoi(argv[++i]); else if (!strcmp(argv[i], "-m")) minlen = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-am")) am = atoi(argv[++i]); else if (!strcmp(argv[i], "-bm")) bm = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-ak")) ak = atoi(argv[++i]); else die("unknown demux argument");
  }
  if (!in || !bcf || !outdir) die("demux needs -i, -b and -o");
  if ((int)strlen(adapter) < ak) ak = strlen(adapter);
  /* barcodes */
  FILE *bf = fopen(bcf, "r"); if (!bf) die("cannot open barcode file");
  sample_t *S = NULL; char **bcs = NULL; int ns = 0, bl = -1; char *line = NULL; size_t lcap = 0; ssize_t ll;
  while ((ll = getline(&line, &lcap, bf)) > 0) {
    while (ll > 0 && (line[ll - 1] == '\n' || line[ll - 1] == '\r')) line[--ll] = 0;
    if (!ll) continue;
    char *tab = strchr(line, '\t'); if (!tab) die("barcode file lines must be <sample>\\t<barcode>");
    *tab = 0; char *bc = tab + 1;
    if (bl < 0) bl = strlen(bc); else if ((int)strlen(bc) != bl) die("all barcodes must have the same length");
    S = xrealloc(S, (ns + 1) * sizeof(sample_t)); bcs = xrealloc(bcs, (ns + 1) * sizeof(char *));
    S[ns].name = strdup(line); bcs[ns] = strdup(bc); S[ns].exact = S[ns].rescued = S[ns].shortn = S[ns].len_sum = S[ns].raw_len_sum = 0; ns++;
  }
  fclose(bf); if (!ns) die("no barcodes");
  smap exact; smap_init(&exact, 1024); for (int s = 0; s < ns; s++) { if (smap_get(&exact, bcs[s], bl) >= 0) die("duplicated barcode"); smap_put(&exact, bcs[s], bl, s); }
  /* 1-mismatch neighbours; ambiguous neighbours (close to two barcodes) map to -2 and are rejected */
  smap near; smap_init(&near, 1 << 16); const char nt[] = "ACGTN"; char tmp[256];
  if (bm >= 1) for (int s = 0; s < ns; s++) for (int p = 0; p < bl; p++) for (int k = 0; k < 5; k++) {
    if (bcs[s][p] == nt[k]) continue;
    memcpy(tmp, bcs[s], bl); tmp[p] = nt[k];
    if (smap_get(&exact, tmp, bl) >= 0) continue;
    int prev = smap_get(&near, tmp, bl); smap_put(&near, tmp, bl, (prev == -1) ? s : -2);
  }
  for (int s = 0; s < ns; s++) {
    char fn[4096]; snprintf(fn, sizeof fn, "%s/%s.%s.seq", outdir, S[s].name, tag);
    S[s].out = fopen(fn, "w"); if (!S[s].out) die("cannot open output file (does outdir exist?)");
    setvbuf(S[s].out, NULL, _IOFBF, 1 << 20);
  }
  /* input */
  FILE *fin; int piped = 0; size_t n = strlen(in);
  if (n > 3 && !strcmp(in + n - 3, ".gz")) {
    char cmd[8192]; snprintf(cmd, sizeof cmd, "(command -v pigz >/dev/null && pigz -dc -p 4 '%s' || gzip -dc '%s')", in, in);
    fin = popen(cmd, "r"); piped = 1;
  } else fin = fopen(in, "r");
  if (!fin) die("cannot open input");
  setvbuf(fin, NULL, _IOFBF, 1 << 22);
  uint64_t total = 0, no_adapter = 0, short_insert = 0, bad_barcode = 0, ambiguous = 0, kept = 0;
  const int minstart = u5 + u3 + bl;   /* earliest adapter position (empty insert); short inserts are counted separately */
  char *h = NULL, *seq = NULL, *plus = NULL, *qual = NULL; size_t hc = 0, sc = 0, pc = 0, qc = 0;
  /* unmatched barcode tally (top list in stats) */
  smap unk; smap_init(&unk, 1 << 20); uint64_t *unkn = calloc(1 << 20, sizeof(uint64_t)); size_t nunk = 0;
  while (getline(&h, &hc, fin) > 0) {
    ssize_t L = getline(&seq, &sc, fin); if (L <= 0) break;
    if (getline(&plus, &pc, fin) <= 0 || getline(&qual, &qc, fin) <= 0) break;
    if (seq[L - 1] == '\n') seq[--L] = 0;
    total++;
    /* scan adapter candidates left to right; the first one preceded by a valid barcode wins, so a chance
       adapter-like k-mer inside the insert does not discard the read */
    int hit = -1, first_hit = -1, s = -1, resc = 0, amb = 0;
    for (int p = minstart; p <= L - 6; p++) {
      int ok;
      if (p <= L - ak) { int mm = 0; for (int j = 0; j < ak && mm <= am; j++) mm += (seq[p + j] != adapter[j]); ok = (mm <= am); }
      else ok = !memcmp(seq + p, adapter, L - p);          /* adapter cut by the read end: exact partial >= 6 nt */
      if (!ok) continue;
      if (first_hit < 0) first_hit = p;
      const char *bc = seq + p - bl;
      int e = smap_get(&exact, bc, bl);
      if (e >= 0) { hit = p; s = e; resc = 0; break; }
      if (bm >= 1) { int q = smap_get(&near, bc, bl); if (q >= 0) { hit = p; s = q; resc = 1; break; } if (q == -2) amb = 1; }
    }
    if (hit < 0) {
      if (first_hit < 0) { no_adapter++; continue; }
      if (amb) { ambiguous++; continue; }
      bad_barcode++;
      const char *bc = seq + first_hit - bl;
      if (nunk < (1 << 19)) { int k = smap_get(&unk, bc, bl); if (k < 0) { smap_put(&unk, bc, bl, nunk); unkn[nunk++] = 1; } else unkn[k]++; }
      continue;
    }
    int ilen = hit - bl - u3 - u5;
    S[s].raw_len_sum += L;
    if (ilen < minlen) { short_insert++; S[s].shortn++; continue; }
    S[s].len_sum += ilen;
    fwrite(seq + u5, 1, ilen, S[s].out); fputc('\n', S[s].out);
    if (resc) S[s].rescued++; else S[s].exact++;
    kept++;
  }
  if (piped) { if (pclose(fin) != 0) die("decompression failed"); } else fclose(fin);
  char fn[4096]; snprintf(fn, sizeof fn, "%s/stats.%s.tsv", outdir, tag);
  FILE *st = fopen(fn, "w"); if (!st) die("cannot write stats");
  /* per sample: reads kept (exact / 1-mismatch barcode), reads with the barcode but a too short insert, summed kept insert
     lengths and summed raw read lengths (all barcode-assigned reads); global rows carry their count in reads_exact */
  fprintf(st, "key\tsample\tbarcode\treads_exact\treads_rescued\treads_short\tinsert_len_sum\traw_len_sum\n");
  for (int s = 0; s < ns; s++) { fclose(S[s].out); fprintf(st, "sample\t%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\n", S[s].name, bcs[s],
    (unsigned long long)S[s].exact, (unsigned long long)S[s].rescued, (unsigned long long)S[s].shortn, (unsigned long long)S[s].len_sum, (unsigned long long)S[s].raw_len_sum); }
  fprintf(st, "total\t\t\t%llu\t0\t0\t0\t0\nno_adapter\t\t\t%llu\t0\t0\t0\t0\nshort_insert\t\t\t%llu\t0\t0\t0\t0\nunknown_barcode\t\t\t%llu\t0\t0\t0\t0\nambiguous_barcode\t\t\t%llu\t0\t0\t0\t0\nkept\t\t\t%llu\t0\t0\t0\t0\n",
          (unsigned long long)total, (unsigned long long)no_adapter, (unsigned long long)short_insert,
          (unsigned long long)bad_barcode, (unsigned long long)ambiguous, (unsigned long long)kept);
  /* top 20 unknown barcodes */
  for (int r = 0; r < 20; r++) {
    size_t best = (size_t)-1; uint64_t bv = 0;
    for (size_t i = 0; i < unk.cap; i++) if (unk.key[i] && unkn[unk.val[i]] > bv) { bv = unkn[unk.val[i]]; best = i; }
    if (best == (size_t)-1) break;
    fprintf(st, "top_unknown_barcode\t\t%s\t%llu\t0\t0\t0\t0\n", unk.key[best], (unsigned long long)bv); unkn[unk.val[best]] = 0;
  }
  fclose(st);
  fprintf(stderr, "demux %s: %llu reads, %llu kept (%.1f%%), no adapter %llu, short %llu, unknown barcode %llu, ambiguous %llu\n", tag,
          (unsigned long long)total, (unsigned long long)kept, total ? 100.0 * kept / total : 0, (unsigned long long)no_adapter,
          (unsigned long long)short_insert, (unsigned long long)bad_barcode, (unsigned long long)ambiguous);
  return 0;
}

/* ---------- collapse ---------- */
typedef struct { char *s; uint32_t len; uint64_t n; } entry;
static entry *T; static size_t tcap, tn;
static void tgrow(void);
static void tadd(const char *s, uint32_t len) {
  if (tn * 2 >= tcap) tgrow();
  size_t i = hash_bytes(s, len) & (tcap - 1);
  while (T[i].s) { if (T[i].len == len && !memcmp(T[i].s, s, len)) { T[i].n++; return; } i = (i + 1) & (tcap - 1); }
  T[i].s = xmalloc(len); memcpy(T[i].s, s, len); T[i].len = len; T[i].n = 1; tn++;
}
static void tgrow(void) {
  size_t oc = tcap; entry *O = T; tcap = oc ? oc * 2 : (1 << 20); T = calloc(tcap, sizeof(entry)); if (!T) die("out of memory");
  for (size_t i = 0; i < oc; i++) if (O[i].s) { size_t j = hash_bytes(O[i].s, O[i].len) & (tcap - 1); while (T[j].s) j = (j + 1) & (tcap - 1); T[j] = O[i]; }
  free(O);
}
static int cmp_count(const void *a, const void *b) { const entry *x = a, *y = b; return (x->n < y->n) - (x->n > y->n); }

static int run_collapse(int argc, char **argv) {
  const char *out = NULL; int minlen = 0, first = -1;
  for (int i = 2; i < argc; i++) {
    if (!strcmp(argv[i], "-o")) out = argv[++i]; else if (!strcmp(argv[i], "-m")) minlen = atoi(argv[++i]);
    else { first = i; break; }
  }
  if (!out || first < 0) die("collapse needs -o and input files");
  char *line = NULL; size_t cap = 0; ssize_t L; uint64_t reads = 0;
  for (int f = first; f < argc; f++) {
    FILE *fi = fopen(argv[f], "r"); if (!fi) die("cannot open collapse input"); setvbuf(fi, NULL, _IOFBF, 1 << 22);
    while ((L = getline(&line, &cap, fi)) > 0) { if (line[L - 1] == '\n') line[--L] = 0; if (L >= minlen && L > 0) { tadd(line, L); reads++; } }
    fclose(fi);
  }
  entry *E = xmalloc((tn ? tn : 1) * sizeof(entry)); size_t k = 0;
  for (size_t i = 0; i < tcap; i++) if (T[i].s) E[k++] = T[i];
  qsort(E, k, sizeof(entry), cmp_count);
  char tmpf[4096]; snprintf(tmpf, sizeof tmpf, "%s.tmp", out);
  gzFile g = gzopen(tmpf, "wb6"); if (!g) die("cannot write collapse output"); gzbuffer(g, 1 << 20);
  for (size_t i = 0; i < k; i++) { gzprintf(g, ">seq%zu_x%llu\n", i + 1, (unsigned long long)E[i].n); gzwrite(g, E[i].s, E[i].len); gzputc(g, '\n'); }
  if (gzclose(g) != Z_OK) die("gzclose failed");
  if (rename(tmpf, out)) die("rename failed");
  fprintf(stderr, "collapse %s: %llu reads, %zu unique\n", out, (unsigned long long)reads, k);
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 2) die("usage: demultiplexer demux|collapse ...");
  if (!strcmp(argv[1], "demux")) return run_demux(argc, argv);
  if (!strcmp(argv[1], "collapse")) return run_collapse(argc, argv);
  die("first argument must be demux or collapse");
  return 2;
}
