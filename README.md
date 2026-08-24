# Applied Genomics and Health AI Projects

A progressive, reproducible portfolio of practical genomics, bioinformatics,
healthcare data-analysis, and artificial-intelligence projects using genuine,
publicly accessible datasets.

## Portfolio objective

Build job-ready skills in biological data acquisition, quality control,
reproducible programming, statistical analysis, visualization, and scientifically
responsible interpretation.

Every project identifies its original data source, states relevant reuse terms,
and documents how its findings were produced. No invented or simulated datasets
are substituted for genuine source data.

## Learning progression

This repository is a progressive curriculum, not a collection of permanently
beginner-level exercises. Explanations will remain beginner-friendly, while the
technical difficulty, independence, and project scope increase over time.

| Phase | Approximate projects | Level | Skills and project direction |
| --- | --- | --- | --- |
| 1. Foundations | 01–04 | Beginner | Python and R basics, biological file formats, public-data access, quality checks, summaries, and plots |
| 2. Applied analysis | 05–09 | Beginner-to-intermediate | Genomic annotation, expression data, statistical testing, reusable functions, and reproducible reports |
| 3. Analytical pipelines | 10–14 | Intermediate-to-advanced | RNA-seq or GWAS workflows, dimensionality reduction, predictive modeling, validation, and explainability |
| 4. Job-ready capstones | 15 onward | Advanced | End-to-end genomics and healthcare AI projects with workflow automation, testing, documentation, ethics, and stakeholder communication |

Progression is competency-based: later exercises will revisit foundational ideas
briefly when needed, but each new project will introduce more realistic data,
stronger analytical reasoning, and more professional engineering practices.

## Completed exercises

| No. | Exercise | Language | Public data source | Main skills |
| --- | --- | --- | --- | --- |
| 01 | [GC-content analysis of cardiovascular and folate-pathway transcripts](genomics-practice-01-python-gc-content-refseq/) | Python | NCBI RefSeq: MTHFR, ACE, AGT, and NOS3 | API access, FASTA parsing, sequence quality control, nucleotide counting, CSV reporting, and scientific interpretation |

## Repository organization

Each exercise is stored in its own numbered project folder:

```text
genomics-practice-01-python-gc-content-refseq/
```

An exercise contains its executable analysis, an explanation of the public data,
reproduction instructions, and small real-data outputs suitable for review.

## Getting started

Install Python 3.10 or newer. Exercise 01 uses only the Python standard library,
so no additional packages are required.

```bash
cd genomics-practice-01-python-gc-content-refseq
python analyze_gc_content.py --email your-real-email@example.com
python -m unittest discover -s tests -v
```

Replace the example email with your own valid email address. NCBI recommends
providing a contact email with automated E-utilities requests.

## Data ethics and reproducibility

- Use public reference records rather than identifiable patient data.
- Link directly to the authoritative source for every dataset.
- Prefer versioned accessions so analyses can be reproduced.
- Check whether a reference record is reviewed, current, or suppressed.
- Distinguish exploratory findings from validated clinical conclusions.
- Follow each data provider's access, attribution, and reuse requirements.
