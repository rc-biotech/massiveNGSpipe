# RNA strand normalization

RNA OFST creation now validates biological read orientation from the final BAM for every RNA sample, including single-end archives and paired archives whose selected mate was aligned as single-end. The original BAM is preserved. The pre-alignment paired mate-selection cache remains separate.

Cached decisions live in `aligned/.rna_strand/<Run>_rna_strand.rds`. Normalized OFST provenance lives in `aligned/.rna_strand/ofst_all/` or `ofst_unique/`, with the same `<Run>_rna_strand.rds` naming. Nothing is added alongside BAM/OFST filenames. Once all sample checks finish, `aligned/QC_STATS/rna_strand.csv` reports every sample in that experiment (non-RNA rows, if present, are marked not applicable).

The detector streams 100,000 primary aligned records at a time, uses NH=1 reads, excludes mitochondrial and overlapping opposite-strand exon regions, respects splice blocks and gives each informative read one vote. It reads another chunk if evidence is insufficient, up to a default cap of 1,000,000 records. At least 1,000 informative reads and 80% agreement are required for same/opposite orientation; 40–60% sense is classified as unstranded. Other cases remain unknown and block RNA conversion, with the study table and cached evidence available for diagnosis. These are conservative initial thresholds, not calibrated confidence probabilities. A coordinate-sorted prefix is a bounded diagnostic sample, not a random genome-wide sample.

Increase the cap when the initial scan is unlucky:

```r
config$rna_strand <- list(max_records = 2000000L)
```

Cache reuse checks BAM size/mtime, annotation signature, detector version and settings. It performs no full-BAM hashing or realignment. Do not use the study's cached paired `strandMode` value as a final-BAM flip: the selected mate may already be sense-oriented.

ORFik must include `readBam(read_strand=...)` and the new RNA provenance helpers. The pipeline checks this before conversion. `read_strand="auto"` uses a valid cached BAM decision; `alignment` preserves the physical alignment, `opposite` reverses biological strand, and `unstranded` retains physical strand but creates total unstranded coverage/counts. Paired `strandMode` retains its existing interpretation. OFST import never repeats a reversal. RNA CIGAR collapse refreshes product provenance after rewriting the file.

Both all/unique mapper conversions share a decision. Changing a decision/product invalidates affected OFST, CIGAR-collapse, merged libraries, coverage, bigWig and count stages. Existing sibling sample markers are preserved; affected derived products and count matrices are copied under `.rna_strand/backups/` before replacement. Completed historical studies are not automatically migrated merely by installing the package. Deliberately reopen their OFST stage after reviewing a repair cohort; do not edit flags for active jobs. Any biological atlas/TE refresh remains a separate downstream operation.

Read-only real-BAM pilots found expected reverse orientation for PRJNA707431, PRJNA287879 and PRJNA748710, forward orientation for PRJNA285961, and insufficient evidence for the poorly mapped PRJNA528026 input. Detection took approximately 0.18–0.78 seconds per sample, annotation preparation 0.82 seconds once per reference, and a vectorized strand reversal of 100,000 alignments 0.005 seconds. An isolated synthetic pipeline exercises conversion, both mapper modes, collapse, coverage, cached reuse, stale-input rejection, tables and sibling-preserving invalidation. No existing production BAM/count/atlas products were repaired during implementation.
