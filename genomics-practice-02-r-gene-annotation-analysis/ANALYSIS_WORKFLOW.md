# Analysis Workflow: Comparative Human Gene Annotation and Transcript Structure

## Analytical objective

Integrate genuine human NCBI Gene annotations with previously generated NCBI
RefSeq transcript measurements to compare chromosome position, strand
orientation, exon structure, genomic interval length, and selected transcript
length for MTHFR, ACE, AGT, and NOS3.

## Input datasets

### NCBI Gene annotations

| Gene | NCBI Gene ID | GRCh38 chromosome accession | Source record |
| --- | ---: | --- | --- |
| MTHFR | 4524 | NC_000001.11 | https://www.ncbi.nlm.nih.gov/gene/4524 |
| ACE | 1636 | NC_000017.11 | https://www.ncbi.nlm.nih.gov/gene/1636 |
| AGT | 183 | NC_000001.11 | https://www.ncbi.nlm.nih.gov/gene/183 |
| NOS3 | 4846 | NC_000007.14 | https://www.ncbi.nlm.nih.gov/gene/4846 |

The reference assembly is GRCh38.p14 (`GCF_000001405.40`). A single NCBI
ESummary request retrieves the four genuine Gene records in JSON format:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=4524,1636,183,4846&retmode=json
```

The documented offline snapshot in
`tests/fixtures/verified_ncbi_gene_annotations.csv` contains real values copied
from those authoritative records and preserves an individual source link for
every gene.

### Project 01 RefSeq transcript measurements

The secondary input is the real-data output created by Project 01:

```text
../genomics-practice-01-python-gc-content-refseq/outputs/gc_content_results.csv
```

Required fields are `gene`, `accession`, `length_nt`, and `gc_percent`.

## Record validation

Before analysis, each NCBI Gene result is checked against:

1. Its expected NCBI Gene ID.
2. Its expected official gene symbol.
3. The organism `Homo sapiens`.
4. Its versioned primary GRCh38 chromosome accession.
5. A positive annotated exon count.

The integration step additionally rejects missing or duplicated Project 01 gene
records and nonpositive RefSeq transcript lengths.

## Coordinate normalization

NCBI Gene ESummary reports `ChrStart` and `ChrStop` as zero-based chromosome
positions. Their order reflects biological orientation, so a minus-strand gene
has `ChrStart > ChrStop`.

The normalization method is:

```text
genomic_start_1_based = min(ChrStart, ChrStop) + 1
genomic_end_1_based   = max(ChrStart, ChrStop) + 1
genomic_span_bp       = abs(ChrStop - ChrStart) + 1
strand                = "+" if ChrStart <= ChrStop; otherwise "-"
```

For example, the genuine MTHFR Gene record has the zero-based biological
orientation:

```text
ChrStart = 11,805,963
ChrStop  = 11,785,722
```

The resulting standardized coordinates are:

```text
Chromosome 1: 11,785,723-11,805,964
Strand: -
Genomic span: 20,242 bp
```

The genuine ACE record illustrates the plus-strand case:

```text
ChrStart = 63,477,060
ChrStop  = 63,498,372

Chromosome 17: 63,477,061-63,498,373
Strand: +
Genomic span: 21,313 bp
```

## Cross-dataset integration

NCBI Gene annotations and Project 01 RefSeq results are joined by official gene
symbol. Original record identifiers remain in the final output:

- NCBI Gene ID identifies the gene-level annotation.
- Versioned chromosome accession identifies the genomic reference sequence.
- Ensembl gene ID provides an independent public cross-reference.
- Versioned RefSeq accession identifies the selected transcript sequence.

The resulting table preserves both gene-level genomic features and transcript-
level measurements, making the distinction between their biological scopes
explicit.

## Derived structural metrics

### Genomic-to-transcript length ratio

```text
genomic_to_transcript_ratio = genomic_span_bp / transcript_length_nt
```

For AGT:

```text
43,061 / 2,148 = 20.05
```

This is a descriptive comparison between the gene's complete annotated genomic
interval and one selected processed RNA reference.

### Transcript fraction of genomic span

```text
transcript_fraction_percent = transcript_length_nt / genomic_span_bp * 100
```

For AGT:

```text
2,148 / 43,061 * 100 = 4.99%
```

This value is not an experimentally measured exon-coverage percentage.

### Annotated exon density

```text
exons_per_10kb = exon_count / genomic_span_bp * 10,000
```

This metric summarizes the reported gene-level exon count relative to the
annotated genomic span; it is not a transcript-specific splice model.

## Visualization

A grouped SVG bar chart compares the GRCh38 genomic interval and selected
RefSeq transcript length for each gene. Both values are expressed in nucleotide
bases, but they represent different biological entities and are interpreted
accordingly.

## Verification strategy

The base-R automated checks verify:

1. The four authentic, individually attributable NCBI Gene records.
2. Real MTHFR minus-strand coordinate normalization.
3. Real ACE plus-strand coordinate normalization.
4. Real AGT and NOS3 inclusive genomic spans.
5. Rejection of invalid chromosome coordinates.
6. Inclusion of all four actual Gene IDs in the NCBI API request.
7. Integration with the genuine Project 01 RefSeq measurements.
8. Detection of missing required Project 01 genes.

No artificial nucleotide sequence or invented genomic annotation is used as an
analysis dataset or test fixture.

## Execution

```bash
Rscript analyze_gene_annotations.R --email your-real-email@example.com
Rscript tests/test_gene_annotations.R
```

To reproduce the analysis using the explicitly documented genuine NCBI Gene
snapshot without requiring a network request or the `jsonlite` package:

```bash
Rscript analyze_gene_annotations.R --offline-fixture
```

## Interpretation boundaries

Gene-level exon counts do not necessarily describe the selected transcript
isoform. Genomic-to-transcript ratios do not quantify expression, intron count,
disease association, or clinical risk. Results are specific to the identified
reference assembly, chromosome accessions, Gene annotations, and RefSeq
transcript versions.
