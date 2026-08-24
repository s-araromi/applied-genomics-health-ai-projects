# Project 01: GC-Content Analysis of Human Cardiovascular and Folate-Pathway Transcripts

## Project overview

This project examines the nucleotide composition and GC content of four human
reference transcripts involved in blood-pressure regulation, vascular function,
and folate metabolism.

The analysis retrieves authentic, versioned transcript sequences from the
National Center for Biotechnology Information Reference Sequence database,
performs sequence quality-control checks, calculates nucleotide frequencies and
GC percentages, and exports reproducible results.

## Analytical question

How do the reference transcripts of MTHFR, ACE, AGT, and NOS3 compare in their
nucleotide composition and GC content?

## Biological context

GC content describes the proportion of guanine and cytosine nucleotides within
a biological sequence. Variation in GC content can influence sequence properties
relevant to molecular biology workflows, including preliminary considerations
for amplification and sequence analysis.

| Gene | Biological role |
| --- | --- |
| MTHFR | Folate metabolism and homocysteine processing |
| ACE | Angiotensin processing and blood-pressure regulation |
| AGT | Angiotensinogen production and blood-pressure regulation |
| NOS3 | Nitric-oxide production and vascular function |

The analysis evaluates reference sequence composition only. It does not measure
gene expression, diagnose disease, establish pathogenicity, or estimate
individual cardiovascular risk.

## Data source

All sequences were obtained from the
[NCBI Reference Sequence database](https://www.ncbi.nlm.nih.gov/refseq/).

| Gene | RefSeq accession | Sequence record |
| --- | --- | --- |
| MTHFR | NM_005957.5 | [View NCBI record](https://www.ncbi.nlm.nih.gov/nuccore/NM_005957.5) |
| ACE | NM_000789.4 | [View NCBI record](https://www.ncbi.nlm.nih.gov/nuccore/NM_000789.4) |
| AGT | NM_001382817.3 | [View NCBI record](https://www.ncbi.nlm.nih.gov/nuccore/NM_001382817.3) |
| NOS3 | NM_000603.5 | [View NCBI record](https://www.ncbi.nlm.nih.gov/nuccore/NM_000603.5) |

The accession suffix identifies the specific version of each sequence and
supports reproducibility.

[Download the four RefSeq transcripts from NCBI EFetch](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_005957.5,NM_000789.4,NM_001382817.3,NM_000603.5&rettype=fasta&retmode=text)

### Reference quality assessment

The older AGT accession `NM_000029.4` is listed as suppressed in the
[NCBI AGT Gene record](https://www.ncbi.nlm.nih.gov/gene/183).

This project uses the reviewed AGT reference transcript `NM_001382817.3`.

### Data usage

NCBI does not impose its own restrictions on the use or redistribution of
molecular database content. Original submitters or countries of origin may
assert additional rights.

See the [NCBI data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
for complete terms.

## Methods

1. Identify and verify four versioned human RefSeq transcript accessions.
2. Retrieve all sequences through one NCBI EFetch API request.
3. Parse FASTA headers and reconstruct complete transcript sequences.
4. Confirm that every requested accession was returned.
5. Count adenine, cytosine, guanine, and thymine nucleotides.
6. Identify and report ambiguous nucleotide characters.
7. Calculate transcript-specific GC percentages.
8. Apply descriptive GC-content categories.
9. Compare transcript-level and pooled GC-content estimates.
10. Export the original sequences, tabulated results, and a Markdown report.
11. Validate the analysis using automated tests and an authentic MTHFR sequence
    subset.

GC content was calculated as:

```text
GC content (%) = [(G count + C count) / (A + C + G + T)] × 100
```

Ambiguous nucleotide characters are excluded from the denominator and reported
separately.

| GC content | Category |
| --- | --- |
| Below 40% | Low GC |
| 40% to 60% | Moderate GC |
| Above 60% | High GC |

These categories are used for descriptive comparison and do not represent
clinical thresholds.

Additional methodological details are provided in
[ANALYSIS_WORKFLOW.md](ANALYSIS_WORKFLOW.md).

## Results

| Gene | RefSeq accession | Length (nt) | A | C | G | T | GC content | Category |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| MTHFR | NM_005957.5 | 7,018 | 1,498 | 2,003 | 2,005 | 1,512 | 57.11% | Moderate GC |
| ACE | NM_000789.4 | 4,962 | 1,041 | 1,593 | 1,404 | 924 | 60.40% | High GC |
| AGT | NM_001382817.3 | 2,148 | 484 | 608 | 565 | 491 | 54.61% | Moderate GC |
| NOS3 | NM_000603.5 | 4,366 | 862 | 1,415 | 1,317 | 772 | 62.57% | High GC |

- Highest transcript-level GC content: NOS3, 62.57%.
- Lowest transcript-level GC content: AGT, 54.61%.
- Mean transcript-level GC content: 58.67%.
- Pooled, length-weighted GC content: 58.99%.
- Ambiguous nucleotide characters across all transcripts: zero.

## Interpretation

NOS3 had the highest GC content among the four reference transcripts, while
AGT had the lowest. ACE and NOS3 exceeded the descriptive 60% threshold, whereas
MTHFR and AGT fell within the 40%–60% range.

The difference between the mean transcript-level GC content and pooled GC
content reflects differences in transcript length. The transcript-level mean
gives each gene equal weight, while the pooled estimate gives greater weight to
longer sequences.

These findings describe whole-transcript nucleotide composition. They do not
demonstrate differences in disease risk, transcriptional activity, clinical
significance, or primer suitability.

## Output files

| File | Description |
| --- | --- |
| `outputs/refseq_transcripts.fasta` | Original transcript sequences retrieved from NCBI |
| `outputs/gc_content_results.csv` | Nucleotide counts, sequence lengths, ambiguous-base counts, GC percentages, and categories |
| `outputs/gc_content_report.md` | Ranked GC-content results and interpretation notes |

## Reproducibility

### Requirements

- Python 3.10 or newer.
- An internet connection for downloading the NCBI records.
- No external Python packages.

### Run the analysis

```bash
python analyze_gc_content.py --email your-real-email@example.com
```

Alternatively:

```bash
python3 analyze_gc_content.py --email your-real-email@example.com
```

Replace the example address with a valid email address. NCBI recommends
including contact information with automated E-utilities requests.

### Run automated tests

```bash
python -m unittest discover -s tests -v
```

The test suite verifies:

- GC-content calculations using a genuine 120-nucleotide MTHFR sequence.
- Reconstruction of multiline FASTA records.
- Correct handling of GC-category boundaries.
- Rejection of empty sequences.

## Technical competencies demonstrated

- Public genomic data acquisition through the NCBI EFetch API.
- Verification of versioned RefSeq accessions.
- Identification of suppressed database records.
- FASTA parsing and sequence reconstruction.
- Sequence completeness and nucleotide quality checks.
- Python functions, dictionaries, loops, and conditional logic.
- Descriptive biological sequence analysis.
- CSV and Markdown report generation.
- Automated unit testing.
- Reproducible and scientifically responsible interpretation.

## Limitations

This analysis includes four selected reference transcripts and does not
represent all cardiovascular, vascular, or folate-pathway genes.

Whole-transcript GC content does not capture local sequence variation across
specific exons, primer-binding regions, or other smaller genomic intervals.

The analysis does not incorporate patient genotypes, variant annotations,
gene-expression measurements, experimental amplification results, or clinical
outcomes.

The GC-content categories are descriptive labels and should not be treated as
validated diagnostic or experimental thresholds.

## References

- [NCBI Reference Sequence database](https://www.ncbi.nlm.nih.gov/refseq/)
- [NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [NCBI data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [NCBI AGT Gene record](https://www.ncbi.nlm.nih.gov/gene/183)
