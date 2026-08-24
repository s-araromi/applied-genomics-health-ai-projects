# Exercise 01: GC-content analysis of real human RefSeq transcripts

## Topic

Use beginner-friendly Python to download and analyze four genuine human reference
transcripts associated with blood-pressure regulation, vascular function, and
folate metabolism.

GC content is the percentage of known sequence letters that are guanine (`G`) or
cytosine (`C`). The analysis calculates nucleotide composition, flags ambiguous
bases, compares transcript-level GC content, and exports reproducible results.

## Dataset

Source: [NCBI Nucleotide / RefSeq](https://www.ncbi.nlm.nih.gov/refseq/).

| Gene | Biological relevance | Versioned RefSeq record |
| --- | --- | --- |
| MTHFR | Folate and homocysteine metabolism | [NM_005957.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_005957.5) |
| ACE | Angiotensin processing and blood-pressure regulation | [NM_000789.4](https://www.ncbi.nlm.nih.gov/nuccore/NM_000789.4) |
| AGT | Angiotensinogen production and blood-pressure regulation | [NM_001382817.3](https://www.ncbi.nlm.nih.gov/nuccore/NM_001382817.3) |
| NOS3 | Nitric-oxide production and vascular function | [NM_000603.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_000603.5) |

The script downloads all four complete transcript sequences directly from NCBI
in one request. These are authentic public reference sequences, not simulated
data and not identifiable patient records.

Direct FASTA access:

https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_005957.5,NM_000789.4,NM_001382817.3,NM_000603.5&rettype=fasta&retmode=text

### Data-quality note

An older AGT accession, `NM_000029.4`, appears in some online references, but
NCBI currently lists it as suppressed. This project deliberately uses the
reviewed AGT accession `NM_001382817.3` instead:

https://www.ncbi.nlm.nih.gov/gene/183

### Data usage

NCBI places no restrictions of its own on the use or redistribution of molecular
database content, while noting that original submitters or countries of origin
may assert additional rights. Acknowledge NCBI and review its complete policy:

https://www.ncbi.nlm.nih.gov/home/about/policies/

## Research question

How do four cardiovascular- and folate-pathway reference transcripts compare in
their nucleotide composition and overall GC content?

## Approach

1. Define four verified, versioned NCBI RefSeq accessions.
2. Retrieve all sequences in FASTA format with one NCBI EFetch API request.
3. Parse the FASTA headers and sequence lines using Python's standard library.
4. Confirm that every requested accession was actually returned.
5. Count `A`, `C`, `G`, and `T` bases and report any ambiguous letters.
6. Calculate GC percentage using only the four known nucleotide bases.
7. Apply descriptive low, moderate, and high GC categories.
8. Compare the unweighted mean GC content with the pooled, length-weighted value.
9. Save the original FASTA, a CSV results table, and a Markdown summary.
10. Test the calculation against a documented real 120-nucleotide MTHFR subset.

## Requirements

- Python 3.10 or newer.
- An internet connection for the initial NCBI download.
- No external Python packages.

## Run the analysis

From this exercise directory:

```bash
python analyze_gc_content.py --email your-real-email@example.com
```

Replace the example address with a real email address. NCBI recommends including
contact information with automated E-utilities requests. The email is sent to
NCBI with the request but is not written into the output files.

On some computers the command is `python3` instead of `python`:

```bash
python3 analyze_gc_content.py --email your-real-email@example.com
```

Run the automated tests:

```bash
python -m unittest discover -s tests -v
```

## Output files

| File | Purpose |
| --- | --- |
| `outputs/refseq_transcripts.fasta` | Downloaded real NCBI transcript sequences and source headers |
| `outputs/gc_content_results.csv` | Gene-level nucleotide counts, ambiguous-base counts, GC percentages, and categories |
| `outputs/gc_content_report.md` | Ranked summary, highest/lowest GC transcripts, and interpretation cautions |

## Scientific interpretation

GC content describes sequence composition; it does not measure gene expression,
determine whether a variant is pathogenic, or diagnose hypertension.

The 40% and 60% category thresholds are illustrative labels for this educational
comparison. They are not universal clinical or experimental cutoffs. Whole-
transcript GC content can inform preliminary sequence assessment, but actual PCR
primer suitability must be evaluated using the specific candidate primer or
amplicon region, along with melting temperature, specificity, secondary
structure, and laboratory conditions.

## Verified results

| Gene | RefSeq accession | Length (nt) | GC (%) | Category |
| --- | --- | ---: | ---: | --- |
| MTHFR | NM_005957.5 | 7,018 | 57.11 | Moderate GC |
| ACE | NM_000789.4 | 4,962 | 60.40 | High GC |
| AGT | NM_001382817.3 | 2,148 | 54.61 | Moderate GC |
| NOS3 | NM_000603.5 | 4,366 | 62.57 | High GC |

- Highest transcript-level GC content: NOS3, 62.57%.
- Lowest transcript-level GC content: AGT, 54.61%.
- Unweighted mean transcript-level GC content: 58.67%.
- Pooled, length-weighted GC content: 58.99%.
- Ambiguous bases identified across all four sequences: zero.

See [WALKTHROUGH.md](WALKTHROUGH.md) for the complete beginner-focused,
step-by-step explanation.

## What I learned

- Find and verify public, versioned NCBI RefSeq records.
- Recognize and avoid a suppressed reference accession.
- Retrieve real biological data through an HTTP API.
- Interpret FASTA headers and reconstruct multi-line sequences.
- Use Python dictionaries, loops, functions, conditionals, and string methods.
- Perform simple sequence quality-control checks.
- Calculate and compare nucleotide-composition statistics.
- Export reproducible CSV and Markdown results.
- Explain biological findings without overstating clinical significance.
- Verify analysis logic with standard-library unit tests.

## References

- NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
- NCBI E-utilities documentation: https://www.ncbi.nlm.nih.gov/books/NBK25499/
- NCBI data-use policies: https://www.ncbi.nlm.nih.gov/home/about/policies/
- NCBI AGT Gene record: https://www.ncbi.nlm.nih.gov/gene/183
