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

## Integration (done, 2026-10-28)
- `pipeline_config(fast_trim = TRUE)` (default): single-end runs with an explicit adapter (or "disable") are trimmed and
  collapsed in one pass by `trimcollapse_run()` (R/trimcollapse.R); outputs keep the pipeline's names and places:
  `trim/report_<Run>.json`, `trim/SINGLE/collapsed_trimmed_<Run>.fasta.gz`; the per-sample `collapsed` marker is set
  and `orient_collapsed_fasta()` applied, so pipeline_collapse() skips these runs. No trimmed fastq is written.
- Covered paths: plain adapter trimming (detect_adapter_and_trim), barcode/UMI trimming with front/tail cuts and the
  post-adapter tail case (fastp_post_adapter_tail semantics, incl. its stdin stage having polyG off), the 1 M-read
  barcode-detection sample (quality filter on), adapter disabled. Fallback to fastp: paired-end, fasta input, fastp
  "auto" adapter detection, polyN special case, anchored-front trimming, missing input.
- The binary compiles on first use; ISA-L when libisal is installed, else zlib (still ~6x).

## Validation
- `tests/testthat/test-trimcollapse.R`: names/locations, ORFik::trimming.table on the report, enable rules, exact
  agreement with fastp 1.3.7 for main pass, post-adapter tail and detection sample (synthetic reads incl. polyG,
  no-adapter and adapter-like-start reads).
- Test matrix (2 M reads x 3 libraries: NextSeq untrimmed, NovaSeq provider-trimmed, UMI+barcode library): main pass,
  front/tail cuts, adapter disabled, detection pass: identical reads, filtering, summaries, adapter counts, curves;
  post-adapter tail identical to fastp_post_adapter_tail(); plain / gz / multi-member gz / zlib build identical;
  truncated gzip rejected.
- `inst/trimcollapse/validate_pipeline_trim.R`: runs the real pipe_trim_collapse() with fast_trim FALSE and TRUE on the
  given fastq files and compares collapsed reads, reports, adapter_barcode_table, markers and flags. Rerun it on a new
  server: `Rscript inst/trimcollapse/validate_pipeline_trim.R <massiveNGSpipe dir> sample1.fastq.gz sample2.fastq.gz`.
- Next speedups: replace fastqc adapter detection (25 s/sample) by counting the candidate adapters inside the tool;
  ISA-L reader for the demultiplexer.
