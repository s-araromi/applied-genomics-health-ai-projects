# Analysis Workflow: Human RefSeq Transcript GC-Content Assessment

## 1. Analysis objective

Assess and compare the nucleotide composition of four human reference
transcripts associated with folate metabolism, blood-pressure regulation, and
vascular function.

The analysis focuses on MTHFR, ACE, AGT, and NOS3 and uses versioned reference
sequences retrieved directly from NCBI RefSeq.

## 2. Input data

| Gene | RefSeq accession | Biological relevance |
| --- | --- | --- |
| MTHFR | NM_005957.5 | Folate metabolism and homocysteine processing |
| ACE | NM_000789.4 | Angiotensin processing and blood-pressure regulation |
| AGT | NM_001382817.3 | Angiotensinogen production |
| NOS3 | NM_000603.5 | Nitric-oxide production and vascular function |

Each accession includes a version suffix identifying the exact sequence used.

- [MTHFR: NM_005957.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_005957.5)
- [ACE: NM_000789.4](https://www.ncbi.nlm.nih.gov/nuccore/NM_000789.4)
- [AGT: NM_001382817.3](https://www.ncbi.nlm.nih.gov/nuccore/NM_001382817.3)
- [NOS3: NM_000603.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_000603.5)

The previous AGT accession `NM_000029.4` was excluded because it is listed as
suppressed in the
[NCBI AGT Gene record](https://www.ncbi.nlm.nih.gov/gene/183).

## 3. Computational environment

The analysis was implemented in Python using the standard library.

| Module | Purpose |
| --- | --- |
| `argparse` | Process command-line arguments |
| `csv` | Export tabulated analytical results |
| `pathlib` | Manage project and output file paths |
| `statistics` | Calculate mean transcript-level GC content |
| `urllib` | Construct and submit the NCBI EFetch request |
| `unittest` | Validate analysis functions |

No third-party packages are required.

## 4. Sequence retrieval

Gene symbols were mapped to versioned RefSeq accessions:

```python
GENE_ACCESSIONS = {
    "MTHFR": "NM_005957.5",
    "ACE": "NM_000789.4",
    "AGT": "NM_001382817.3",
    "NOS3": "NM_000603.5",
}
```

The NCBI EFetch request used these parameters:

```python
parameters = {
    "db": "nuccore",
    "id": ",".join(GENE_ACCESSIONS.values()),
    "rettype": "fasta",
    "retmode": "text",
    "tool": "applied_genomics_health_ai_projects",
}
```

| Parameter | Value | Purpose |
| --- | --- | --- |
| `db` | `nuccore` | Select the NCBI nucleotide database |
| `id` | Four comma-separated accessions | Identify the requested records |
| `rettype` | `fasta` | Request FASTA-formatted sequences |
| `retmode` | `text` | Return records as plain text |
| `tool` | `applied_genomics_health_ai_projects` | Identify the requesting application |

[Retrieve the four transcript sequences](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_005957.5,NM_000789.4,NM_001382817.3,NM_000603.5&rettype=fasta&retmode=text)

## 5. FASTA parsing

The parser:

1. Reads the downloaded FASTA content line by line.
2. Detects record headers beginning with `>`.
3. Extracts accession identifiers.
4. Combines multiline nucleotide entries into complete sequences.
5. Converts sequence characters to uppercase.
6. Stores each sequence under its corresponding accession.

The parser rejects sequence content encountered before a valid FASTA header.

## 6. Sequence quality assessment

The analysis verifies that all requested accession numbers are present in the
downloaded FASTA data.

If any record is absent, execution stops and identifies the missing accession.

Each reconstructed sequence is checked for recognizable nucleotide characters
before calculating GC content.

## 7. Nucleotide composition

```python
base_counts = {
    base: sequence.count(base)
    for base in "ACGT"
}
```

```python
known_base_count = sum(base_counts.values())
```

```python
ambiguous_base_count = len(sequence) - known_base_count
```

Although the selected records represent messenger RNA transcripts, NCBI FASTA
sequences use the DNA-style nucleotide alphabet `A`, `C`, `G`, and `T`.

Ambiguous characters are reported separately and excluded from the denominator
used for GC-content calculations.

## 8. GC-content calculation

```text
GC content (%) = [(G + C) / (A + C + G + T)] × 100
```

```python
gc_count = base_counts["G"] + base_counts["C"]

gc_percent = round(
    (gc_count / known_base_count) * 100,
    2,
)
```

| GC percentage | Classification |
| --- | --- |
| Less than 40% | Low GC |
| 40% through 60% | Moderate GC |
| Greater than 60% | High GC |

The categories support descriptive comparison and do not represent clinical
decision thresholds.

## 9. Summary statistics

### Mean transcript-level GC content

Each transcript contributes equally.

```text
Mean transcript-level GC content = 58.67%
```

### Pooled GC content

All recognized nucleotide bases are combined before calculating the overall
percentage.

```text
Pooled GC content = 58.99%
```

The pooled estimate gives greater weight to longer transcript sequences.

## 10. Output files

| File | Description |
| --- | --- |
| `outputs/refseq_transcripts.fasta` | Original transcript sequence records |
| `outputs/gc_content_results.csv` | Nucleotide composition and GC-content results |
| `outputs/gc_content_report.md` | Ranked findings and interpretation notes |

## 11. Analytical results

| Gene | RefSeq accession | Length (nt) | GC content | Classification |
| --- | --- | ---: | ---: | --- |
| MTHFR | NM_005957.5 | 7,018 | 57.11% | Moderate GC |
| ACE | NM_000789.4 | 4,962 | 60.40% | High GC |
| AGT | NM_001382817.3 | 2,148 | 54.61% | Moderate GC |
| NOS3 | NM_000603.5 | 4,366 | 62.57% | High GC |

```text
Highest GC content: NOS3, 62.57%.
Lowest GC content: AGT, 54.61%.
Mean transcript-level GC content: 58.67%.
Pooled GC content: 58.99%.
Ambiguous bases: 0.
```

## 12. Automated validation

The automated tests are located in:

```text
tests/test_gc_content.py
```

The test suite uses the first 120 genuine nucleotide bases from the MTHFR
transcript `NM_005957.5`.

[Retrieve the first 120 bases of MTHFR NM_005957.5](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_005957.5&seq_start=1&seq_stop=120&rettype=fasta&retmode=text)

| Nucleotide | Count |
| --- | ---: |
| A | 25 |
| C | 39 |
| G | 42 |
| T | 14 |
| Total | 120 |

```text
GC content = [(39 + 42) / 120] × 100

GC content = 67.50%
```

Additional tests validate multiline FASTA reconstruction, classification
boundaries, and rejection of empty sequences.

```bash
python -m unittest discover -s tests -v
```

## 13. Interpretation and limitations

NOS3 has the highest GC content among the selected transcripts, while AGT has
the lowest. ACE and NOS3 exceed the descriptive 60% threshold.

Reference transcript composition does not establish disease risk, patient
genotype, gene-expression levels, or PCR primer suitability.

The four selected genes do not represent all cardiovascular or folate-related
genomic pathways.

## 14. Reproduction

```bash
cd genomics-practice-01-python-gc-content-refseq
```

```bash
python analyze_gc_content.py --email your-real-email@example.com
```

```bash
python -m unittest discover -s tests -v
```

## References

- [NCBI Reference Sequence database](https://www.ncbi.nlm.nih.gov/refseq/)
- [NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [NCBI AGT Gene record](https://www.ncbi.nlm.nih.gov/gene/183)
- [NCBI data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
