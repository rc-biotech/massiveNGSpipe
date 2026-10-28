# trimcollapse: fused fastp-compatible trimming + collapsing (branch fast-trim-collapse, 2026-10-28)

## Result
20 M reads (Yewdell Ribo-seq, 201 nt NextSeq 2000 reads, TruSeq adapter in 97%), on the shared server under load
(load average 25-43 of 48 cores, single runs):

| Path | Wall time | Notes |
|---|---|---|
| current: fastqc adapter detection + fastp 1.3.7 (16 thr) + ORFik::collapse.fastq | ~245 s (25 + 86 + 134) | 3.3 GB uncompressed trimmed fastq, 4.6 GB R memory |
| fastp writing to stdout (no file) | 58 s | trimming only, no collapse |
| trimcollapse v1, single thread, zlib | 194 s | |
| trimcollapse v2, zlib reader, 4-16 workers | 38-41 s | decompression bound (zlib) |
| trimcollapse v2, ISA-L reader, sharded collapse, pigz output, 8-24 workers | **15.6-20 s** | 12-16x; no temp file |
| floor: ISA-L decompression alone | 7.7 s | pigz -dc: 41 s, gzip -dc: 48 s |

## Faithfulness (verified)
- Trimmed reads: identical to fastp 1.3.7 for every read (2 M reads, two libraries: NextSeq 2000 untrimmed and NovaSeq
  provider-trimmed), including fastp's one-insertion/one-deletion adapter fallback (which compares the read start for
  every position, a fastp quirk that is reproduced on purpose), polyG auto-trimming for two-colour instruments, the
  quality/N/length filters in fastp's order and adapter-dimer classification (--dimer_max_len 2).
- Collapsed output: identical to fastp + ORFik::collapse.fastq on 20 M reads (7,149,756 unique sequences, same counts).
- JSON report: same structure and keys as fastp 1.3.7; identical summary before/after, filtering_result,
  adapter_trimmed_reads/bases, read1_before/after_filtering totals, quality/content curves and kmer counts.
  Differences: read1_adapter_counts groups the cut adapter tails slightly differently (top 10 + others);
  duplication.rate is exact from the collapse (fastp estimates it from a raw-read sample); no html report;
  fastp_version names this tool; the command line records the real input path (input is read directly, not streamed).
  None of these differing fields is read by massiveNGSpipe or ORFik (they use summary, adapter_cutting sequence/reads,
  total_cycles and the after-filtering curves).

## Design
reader thread (ISA-L igzip inflate, zlib fallback; multi-member gzip; record-aligned 8 MB blocks) -> N workers
(parse, fastp statistics, trim, filter; hash each passing read) -> 8 collapse shards (one hash table per thread, no
locks) -> sort by count -> pigz output (zlib fallback). Output and report do not depend on the worker count.
Single-end / read 1 only (paired-end collapse uses read 1 in the pipeline too; the PE fastp report keeps fastp for now).

## Integration plan (not on master yet)
1. R wrapper `trim_collapse_fast()` (compile on first use like the demultiplexer; -DHAVE_ISAL when libisal is present).
2. In pipeline_trim(): optional `config$fast_trim` (default FALSE at first): barcode-detection pass = trimcollapse
   `-r 1000000 --fastq-out` (same JSON + trimmed sample, so barcode detection code is unchanged); full pass writes
   `trim/report_<Run>.json` and `trim/SINGLE/collapsed_trimmed_<Run>.fasta.gz` directly and sets both the trim and
   collapsed sample markers (no trimmed fastq, no R collapse).
3. Adapter detection: replace fastqc (25 s/sample) by counting the same candidate adapters in the first reads inside
   the tool (< 1 s), validated against fastqc_adapters_info() on the processed studies.
4. Validation before default-on: rerun a set of finished studies (SE Ribo-seq, with/without barcodes and UMIs, several
   instruments) and require identical collapsed files and report summaries.
5. Same ISA-L reader in the demultiplexer (currently pigz-bound: ~0.65 M reads/s per file).
