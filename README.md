# Applied Genomics and Health AI Projects

My progressive, reproducible portfolio of practical genomics, bioinformatics,
healthcare data analysis, and artificial intelligence projects using genuine,
publicly accessible datasets.

## Portfolio objective

Build skills in biological data acquisition, quality control, reproducible
programming, statistical analysis, visualization, and scientifically responsible
interpretation.

Every project identifies its original data source, states relevant reuse terms,
and documents how its findings were produced. No invented or simulated datasets
are substituted for genuine source data.

## Project development

Projects are organized chronologically and increase in analytical and technical
scope. Initial projects focus on biological data retrieval, sequence analysis,
data cleaning, quality assessment, and visualization. Subsequent projects extend
into genomic annotation, gene expression, statistical modelling, RNA sequencing,
genome-wide association studies, machine learning, model evaluation, workflow
automation, and healthcare AI.

Each project documents its analytical question, public dataset, methodology,
reproducible code, results, interpretation, and limitations.

## Projects

| No. | Project | Language | Public data source | Methods and technical skills |
| --- | --- | --- | --- | --- |
| 01 | [GC-content analysis of cardiovascular and folate-pathway transcripts](genomics-practice-01-python-gc-content-refseq/) | Python | NCBI RefSeq: MTHFR, ACE, AGT, and NOS3 | API access, FASTA parsing, sequence quality control, nucleotide counting, CSV reporting, and scientific interpretation |
| 02 | [Comparative genomic annotation and gene-transcript structure](genomics-practice-02-r-gene-annotation-analysis/) | R | NCBI Gene, GRCh38.p14, and Project 01 RefSeq results | Gene annotation, coordinate normalization, strand interpretation, data integration, comparative visualization, and reproducibility checks |
| 03 | [Paired-end FASTQ read quality control](genomics-practice-03-python-fastq-read-quality-control/) | Python | ENA run ERR194147, NA12878 whole-genome sequencing | Streaming FASTQ subsets, paired-read validation, Phred scores, per-cycle quality, checksums, and QC visualization |

## Repository organization

Each project is stored in its own numbered folder:

```text
genomics-practice-01-python-gc-content-refseq/
genomics-practice-02-r-gene-annotation-analysis/
genomics-practice-03-python-fastq-read-quality-control/
```

Each project contains the analysis code, dataset documentation, reproduction
instructions, results, and relevant output files.

## Getting started

Requirements and commands are documented within each project. Project 03 uses
Python 3.10 or newer and `matplotlib`:

```bash
cd genomics-practice-03-python-fastq-read-quality-control
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python analyze_fastq_quality.py
python -m unittest discover -s tests -v
```

## Data ethics and reproducibility

- Use publicly available reference records instead of identifiable patient data.
- Link directly to the authoritative source for every dataset.
- Use versioned accession numbers where available to support reproducibility.
- Check whether reference records are reviewed, current, or suppressed.
- Distinguish exploratory findings from validated clinical conclusions.
- Follow each data provider's access, attribution, and reuse requirements.
