# Beginner walkthrough: analyzing genuine human RefSeq transcripts

## 1. Understand the biological problem

A gene is a DNA region that contributes to producing a functional RNA or protein.
A transcript is an RNA copy produced from a gene. Messenger RNA, or mRNA, carries
instructions that cells can use to build a protein.

This exercise uses four genuine human reference mRNA records:

- `MTHFR` contributes to folate metabolism and homocysteine processing.
- `ACE` participates in the biochemical pathway that regulates blood pressure.
- `AGT` encodes angiotensinogen, an upstream component of the same pathway.
- `NOS3` supports nitric-oxide production and blood-vessel function.

A reference sequence is a curated public sequence used as a reproducible
scientific reference. Here, each reference comes from NCBI RefSeq.

The research question is: how do these four reference transcripts differ in
their nucleotide composition and overall GC content?

## 2. Define the essential sequence terms

A nucleotide is one letter in a biological sequence. NCBI represents these
transcript sequences using four DNA-style letters:

- `A`: adenine.
- `C`: cytosine.
- `G`: guanine.
- `T`: thymine.

Although mature RNA biologically contains uracil, `U`, instead of thymine, `T`,
RefSeq nucleotide FASTA records commonly represent mRNA sequences using the
DNA-style `A/C/G/T` alphabet.

GC content is the percentage of known sequence letters that are `G` or `C`:

```text
GC percentage = (number of G bases + number of C bases)
                -------------------------------------- x 100
                number of A + C + G + T bases
```

An ambiguous base is a letter whose exact identity is uncertain, such as `N`.
We report such letters separately and exclude them from the GC denominator.

## 3. Install and verify Python

Python is the programming language used to retrieve and analyze the sequences.

Open a terminal. On Windows, use PowerShell, Command Prompt, or the terminal
inside Visual Studio Code. Then enter:

```bash
python --version
```

If your computer recognizes `python3` but not `python`, use:

```bash
python3 --version
```

The result should show Python 3.10 or newer. This exercise needs no additional
packages because it uses the Python standard library: tools bundled with Python.

## 4. Understand the repository layout

A repository is a project folder whose history can be tracked using Git.
Git is software for recording changes, and GitHub is an online service that can
host Git repositories.

The portfolio's root directory is:

```text
applied-genomics-health-ai-projects/
```

The first exercise is:

```text
genomics-practice-01-python-gc-content-refseq/
```

Its main Python file is `analyze_gc_content.py`. Its `outputs` directory stores
real downloaded sequences and analysis results. Its `tests` directory checks
that the calculations work correctly.

## 5. Verify the real public dataset

An accession is a unique identifier assigned to a sequence record. The suffix
after the period identifies the exact version.

For example, `NM_005957.5` means version 5 of the MTHFR transcript accession
`NM_005957`.

Open each original record:

- MTHFR: https://www.ncbi.nlm.nih.gov/nuccore/NM_005957.5
- ACE: https://www.ncbi.nlm.nih.gov/nuccore/NM_000789.4
- AGT: https://www.ncbi.nlm.nih.gov/nuccore/NM_001382817.3
- NOS3: https://www.ncbi.nlm.nih.gov/nuccore/NM_000603.5

An older AGT record, `NM_000029.4`, appears in older material but is currently
listed as suppressed in the NCBI AGT Gene record. A suppressed accession has
been withdrawn from the active reference collection. This project therefore
uses the reviewed replacement `NM_001382817.3`.

The NCBI Gene record showing both the reviewed and suppressed accessions is:

https://www.ncbi.nlm.nih.gov/gene/183

## 6. Read the accession dictionary

A dictionary stores labeled values. Here the gene name is the key and its
versioned accession is the value:

```python
GENE_ACCESSIONS = {
    "MTHFR": "NM_005957.5",
    "ACE": "NM_000789.4",
    "AGT": "NM_001382817.3",
    "NOS3": "NM_000603.5",
}
```

This structure keeps the gene-to-record relationship explicit and reproducible.

## 7. Understand the NCBI API request

An API, or application programming interface, is a defined way for one program
to request information from another system.

The script uses NCBI EFetch, an official API for retrieving database records.

```python
parameters = {
    "db": "nuccore",
    "id": ",".join(GENE_ACCESSIONS.values()),
    "rettype": "fasta",
    "retmode": "text",
    "tool": "applied_genomics_health_ai_projects",
}
```

The parameters mean:

- `db="nuccore"`: search the nucleotide sequence database.
- `id`: request the four specific versioned accessions.
- `rettype="fasta"`: ask for biological sequences in FASTA format.
- `retmode="text"`: return ordinary text.
- `tool`: identify the application making the request.

The expression `",".join(...)` combines four accession values into one
comma-separated string. One combined request is kinder to the public API than
four separate requests.

## 8. Define FASTA format

FASTA is a common plain-text format for biological sequences. Each record starts
with a header line beginning with `>` and is followed by sequence letters.

The project obtains the original FASTA text directly from NCBI. You can inspect
the downloaded records in:

```text
outputs/refseq_transcripts.fasta
```

A parser is code that converts a structured text format into information a
program can work with. Our FASTA parser identifies header lines and joins the
following sequence lines into one uninterrupted sequence.

## 9. Understand counting and quality control

Quality control means checking whether data are complete and suitable for the
intended analysis. This script verifies that all four accessions were returned
and checks whether any unexpected sequence letters appear.

```python
base_counts = {base: sequence.count(base) for base in "ACGT"}
known_base_count = sum(base_counts.values())
ambiguous_base_count = len(sequence) - known_base_count
```

`sequence.count(base)` counts how often a letter occurs. `len(sequence)` returns
the total number of letters. `sum(...)` adds the four known-base counts together.

The first line is a dictionary comprehension: a compact way to build a
dictionary by repeating an operation for each value in a collection. It creates
one count for each of `A`, `C`, `G`, and `T`.

## 10. Calculate and classify GC content

```python
gc_count = base_counts["G"] + base_counts["C"]
gc_percent = round((gc_count / known_base_count) * 100, 2)
```

`round(..., 2)` displays the result to two decimal places. The program uses the
following descriptive labels:

- Below 40%: `Low GC`.
- From 40% through 60%: `Moderate GC`.
- Above 60%: `High GC`.

These thresholds are convenient teaching labels, not validated medical or
experimental decision thresholds.

## 11. Run the full analysis

Navigate into the exercise folder:

```bash
cd applied-genomics-health-ai-projects/genomics-practice-01-python-gc-content-refseq
```

Run the program:

```bash
python analyze_gc_content.py --email your-real-email@example.com
```

Replace `your-real-email@example.com` with your real email address. NCBI
recommends contact information for automated requests. The address is supplied
only while running the command; it is not hardcoded into the project or saved
to the generated result files.

`--email` is a command-line argument: a named value supplied after the program
name to configure one execution.

## 12. Compare your output with the verified results

```text
MTHFR | NM_005957.5    | 7,018 nt | GC: 57.11% | Moderate GC
ACE   | NM_000789.4    | 4,962 nt | GC: 60.40% | High GC
AGT   | NM_001382817.3 | 2,148 nt | GC: 54.61% | Moderate GC
NOS3  | NM_000603.5    | 4,366 nt | GC: 62.57% | High GC

Highest GC content: NOS3 (62.57%)
Lowest GC content: AGT (54.61%)
Mean transcript-level GC content: 58.67%
```

The Markdown report also gives a pooled, length-weighted GC percentage of
58.99%. The unweighted mean gives each transcript equal importance. The pooled
value gives longer transcripts more influence because it combines all bases.

## 13. Inspect the output files

`outputs/refseq_transcripts.fasta` contains the real source records.

`outputs/gc_content_results.csv` contains a CSV, or comma-separated values,
table. CSV files can be opened using Excel, R, Python, or spreadsheet software.

`outputs/gc_content_report.md` contains a Markdown report. Markdown is a simple
text formatting language that GitHub renders as headings, tables, and links.

## 14. Run the tests

A unit test checks a small part of a program against an expected result. The
project tests a genuine 120-nucleotide subset of the MTHFR reference sequence.

```bash
python -m unittest discover -s tests -v
```

The documented MTHFR subset contains 25 A, 39 C, 42 G, and 14 T bases. Its GC
content is therefore `(39 + 42) / 120 x 100 = 67.5%`.

The test command also checks FASTA parsing, category boundaries, and rejection
of an empty sequence. All four tests should report `ok`.

## 15. Interpret the findings carefully

NOS3 has the highest whole-transcript GC percentage in this four-record subset,
while AGT has the lowest. ACE and NOS3 exceed the illustrative 60% threshold.

However:

- Whole-transcript GC does not determine hypertension risk.
- Reference sequences are not patient genotypes.
- GC percentage does not reveal gene-expression levels.
- A whole-transcript average cannot establish whether a specific PCR primer is
  suitable; the actual primer or amplicon must be assessed separately.

An amplicon is a DNA segment amplified during polymerase chain reaction, or PCR.

## 16. Commit the project with Git

From the portfolio root:

```bash
git init
git add README.md .gitignore genomics-practice-01-python-gc-content-refseq
git commit -m "Add real-data RefSeq GC-content analysis portfolio project"
```

A commit is a saved snapshot of your project. To publish it, first create an
empty repository named `applied-genomics-health-ai-projects` on GitHub. Then:

```bash
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/applied-genomics-health-ai-projects.git
git push -u origin main
```

Replace `YOUR-USERNAME` with your own GitHub username. A remote is the online
copy of your repository, and `git push` uploads your local commits to it.

## Troubleshooting

If `python` is not recognized, install Python or use `python3`.

If the download fails, check your internet connection and try the direct FASTA
URL in the exercise README. Temporary NCBI service errors can usually be handled
by retrying after a short wait.

If a requested accession is missing, open its NCBI record to determine whether
it has been updated, replaced, or suppressed. Do not silently substitute a
different sequence: document any replacement accession and rerun the analysis.
