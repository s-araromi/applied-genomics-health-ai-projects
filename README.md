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
| 02 | [Comparative genomic annotation and gene-transcript structure analysis](genomics-practice-02-r-gene-annotation-analysis/) | R | NCBI Gene: MTHFR, ACE, AGT, and NOS3 | Gene coordinates, strand orientation, exon summaries, reference-assembly handling, and cross-project integration |
| 03 | [Paired-end FASTQ read quality control](genomics-practice-03-python-fastq-read-quality-control/) | Python | ENA run ERR194147, NA12878 whole-genome sequencing | Streaming FASTQ subsets, paired-read validation, Phred scores, per-cycle quality, checksums, and QC visualization |
| 04 | [Paired GEO expression-matrix quality assessment](genomics-practice-04-r-geo-expression-quality-control/) | R | NCBI GEO GSE19804 and GPL570 | Expression-matrix parsing, paired-sample validation, distributional QC, correlations, review flags, and reproducible graphics |
| 05 | [Paired lung-tissue differential-expression analysis](genomics-practice-05-python-paired-differential-expression/) | Python | NCBI GEO GSE19804 and GPL570 | Paired statistical testing, BH correction, effect-size criteria, annotation integration, sensitivity analysis, and gene-level reporting |
| 06 | [Functional enrichment of paired lung-tumor expression signals](genomics-practice-06-r-functional-enrichment/) | R | GSE19804-derived results, Gene Ontology, and Reactome | Platform-aware over-representation, ontology propagation, evidence sensitivity, FDR control, redundancy reduction, and pathway reporting |
| 07 | [CFTR ClinVar classification-conflict audit](genomics-practice-07-python-clinvar-conflict-audit/) | Python | NCBI ClinVar, 2026-09-28 snapshot | E-utilities retrieval, JSON/XML parsing, versioned-record validation, conflict decomposition, evidence triage, and reproducible reporting |

## Repository organization

Each project is stored in its own numbered folder:

```text
genomics-practice-NN-language-short-kebab-case-topic/
```

Each project contains the analysis code, dataset documentation, reproduction
instructions, results, and relevant output files.

## Getting started

Requirements and commands are documented within each project. For example,
Project 07 requires Python 3.10 or newer and Matplotlib.

```bash
cd genomics-practice-07-python-clinvar-conflict-audit
python3 -m pip install -r requirements.txt
python3 analyze_clinvar_conflicts.py
python3 -m unittest discover -s tests -v
```

## Data ethics and reproducibility

- Use publicly available reference records instead of identifiable patient data.
- Link directly to the authoritative source for every dataset.
- Use versioned accession numbers where available to support reproducibility.
- Check whether reference records are reviewed, current, or suppressed.
- Distinguish exploratory findings from validated clinical conclusions.
- Follow each data provider’s access, attribution, and reuse requirements.
