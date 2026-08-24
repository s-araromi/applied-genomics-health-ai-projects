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

## Repository organization

Each project is stored in its own numbered folder:

```text
genomics-practice-01-python-gc-content-refseq/
```

Each project contains the analysis code, dataset documentation, reproduction
instructions, results, and relevant output files.

## Getting started

Project 01 requires Python 3.10 or newer and uses only the Python standard
library. No additional packages are required.

```bash
cd genomics-practice-01-python-gc-content-refseq
python analyze_gc_content.py --email your-real-email@example.com
python -m unittest discover -s tests -v
```

Replace the example email address with your own valid email address. NCBI
recommends providing contact information with automated E-utilities requests.

## Data ethics and reproducibility

- Use publicly available reference records instead of identifiable patient data.
- Link directly to the authoritative source for every dataset.
- Use versioned accession numbers where available to support reproducibility.
- Check whether reference records are reviewed, current, or suppressed.
- Distinguish exploratory findings from validated clinical conclusions.
- Follow each data provider’s access, attribution, and reuse requirements.
