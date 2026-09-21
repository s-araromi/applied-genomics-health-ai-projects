# Project 03: Paired-End FASTQ Read Quality Control

## Project summary

This project performs reproducible quality control on the first 5,000 paired-end
reads from the public human whole-genome sequencing run
[`ERR194147`](https://www.ebi.ac.uk/ena/browser/view/ERR194147). The run represents
the NA12878 reference sample and was generated on an Illumina HiSeq 2000.

The workflow streams only the requested leading records from ENA's compressed
FASTQ files. It does not download the complete source files, which together
exceed 100 GB.

## Analytical objectives

- Validate four-line FASTQ structure and Phred+33 quality encoding.
- Confirm that R1 and R2 contain matching read-pair identifiers.
- Calculate read counts, base counts, read lengths, mean Phred quality, Q30
  percentages, GC content, and ambiguous-base counts.
- Summarize quality by sequencing cycle.
- Preserve SHA-256 checksums for the exact local subsets.
- Produce tabular, narrative, and graphical QC outputs.

## Data source

| Item | Value |
| --- | --- |
| Repository | European Nucleotide Archive (ENA) |
| Run accession | `ERR194147` |
| Sample | NA12878 |
| Organism | *Homo sapiens* |
| Assay | Whole-genome sequencing |
| Platform | Illumina HiSeq 2000 |
| Subset | First 5,000 paired reads |

Source files:

- [ERR194147 R1 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz)
- [ERR194147 R2 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz)

The source is public. Reuse remains subject to the
[EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use).

## Verified reference results

These results describe the specified 5,000-pair subset:

| Metric | R1 | R2 | Combined |
| --- | ---: | ---: | ---: |
| Reads | 5,000 | 5,000 | 10,000 |
| Bases | 505,000 | 505,000 | 1,010,000 |
| Mean Phred quality | 35.93 | 35.53 | 35.73 |
| Bases reaching Q30 | 92.80% | 91.53% | 92.16% |
| GC content | 49.87% | 49.74% | 49.81% |
| Ambiguous `N` bases | 12 | 1,292 | 1,304 |

R1 met both declared QC criteria. R2 exceeded the descriptive 0.1% threshold
for ambiguous bases and was marked for review. This is a technical flag for the
subset, not evidence that the complete sequencing run is unusable.

## Requirements

- Python 3.10 or newer.
- Internet access to stream the subset from ENA.
- `matplotlib` for PNG and SVG figures.

## Reproduction

Create and activate an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the analysis:

```bash
python analyze_fastq_quality.py
```

The downloaded subsets are cached in `data/` and excluded from Git. Subsequent
runs reuse them. Existing local FASTQ files can be supplied explicitly:

```bash
python analyze_fastq_quality.py --r1 path/to/R1.fastq.gz --r2 path/to/R2.fastq.gz
```

Run the offline automated tests:

```bash
python -m unittest discover -s tests -v
```

## Outputs

| File | Description |
| --- | --- |
| `outputs/fastq_qc_summary.csv` | Read-level QC metrics and review status |
| `outputs/per_cycle_quality.csv` | Per-cycle mean Phred and Q30 percentages |
| `outputs/subset_checksums.csv` | Byte counts and SHA-256 checksums for local subsets |
| `outputs/fastq_qc_report.md` | Narrative summary and interpretation |
| `outputs/fastq_quality_overview.png` | Raster QC visualization |
| `outputs/fastq_quality_overview.svg` | Vector QC visualization |

## Interpretation limits

The analysis covers only the first 5,000 read pairs. It is suitable for learning
FASTQ structure and implementing reproducible QC, but it does not replace a
complete run-level assessment. Sequence quality does not establish genotype,
variant pathogenicity, disease risk, or clinical validity.
