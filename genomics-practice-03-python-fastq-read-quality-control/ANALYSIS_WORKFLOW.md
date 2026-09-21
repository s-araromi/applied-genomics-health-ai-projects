# Analysis Workflow: Paired-End FASTQ Read Quality Control

## 1. Objective

Assess the technical quality of a precisely defined authentic paired-end FASTQ
subset while avoiding a full download of the very large source run.

## 2. Source selection

ENA run `ERR194147` was selected because it is a public human whole-genome
sequencing dataset for the extensively studied NA12878 reference sample. The
analysis uses the first 5,000 records from each mate file.

## 3. Streaming subset acquisition

`urllib.request` opens each HTTPS source. `gzip.GzipFile` decompresses the
response as a stream, and the workflow stops after exactly 5,000 validated
records. Each subset is written as deterministic gzip content and identified by
a locally calculated SHA-256 checksum.

The complete source archives are not downloaded or committed.

## 4. FASTQ validation

Every record must contain:

1. An identifier line beginning with `@`.
2. A nucleotide sequence containing only `A`, `C`, `G`, `T`, or `N`.
3. A separator line beginning with `+`.
4. A quality string with the same length as the sequence.

Quality characters are decoded as Phred+33. R1 and R2 are also required to
contain the same number of records and matching normalized identifiers in the
same order.

## 5. Metrics

For each mate, the workflow calculates:

- Number of reads and bases.
- Minimum and maximum read length.
- Mean Phred quality across all bases.
- Percentage of bases with Phred score at least 30.
- GC percentage among recognized `A`, `C`, `G`, and `T` bases.
- Count and percentage of ambiguous `N` bases.
- Mean quality and Q30 percentage for every sequencing cycle.

The combined row is weighted by observed base counts rather than averaging the
two mate-level percentages.

## 6. Descriptive review criteria

A mate is marked `PASS` when at least 90% of bases reach Q30 and no more than
0.1% are ambiguous `N` calls. Otherwise it is marked `REVIEW`.

These project criteria support transparent comparison. They are not universal
sequencing acceptance thresholds.

## 7. Outputs

The workflow writes machine-readable CSV files, a Markdown report, a PNG figure,
and an SVG figure. Raw subsets remain in the ignored `data/` directory.

## 8. Automated validation

Ten offline unit tests cover FASTQ parsing, base counts, Phred decoding,
identifier normalization, mate matching, malformed records, invalid bases,
combined calculations, and QC status assignment.

## 9. Limitations

- The leading 5,000 pairs are not a random sample of the run.
- Subset metrics cannot establish whole-run suitability.
- The workflow performs quality assessment, not adapter trimming, alignment,
  duplicate marking, variant calling, or clinical interpretation.

## References

- [ENA record ERR194147](https://www.ebi.ac.uk/ena/browser/view/ERR194147)
- [ENA FASTQ file format guidance](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html)
- [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use)
