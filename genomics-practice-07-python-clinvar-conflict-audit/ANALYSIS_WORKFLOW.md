# Project 07 Analysis Workflow

## Analytical objective

Construct a validated, reproducible review queue for high-submission CFTR variants
whose current aggregate ClinVar germline classification is conflicting. The work
separates data retrieval, source validation, evidence summarization, operational
prioritization, and clinical interpretation.

## Data flow

```text
ClinVar ESearch query
        |
        v
Complete Variation ID set (549)
        |
        v
ESummary screening for exact aggregate conflict (407)
        |
        v
Rank by supporting SCV count; select 25
        |
        v
Individual VCV XML retrieval and current-status validation
        |
        v
Structured fields, conflict groups, quality flags, and priority tiers
        |
        v
CSV tables, report, and overview figure
```

## Source query and subset definition

The 2026-09-28 snapshot used this exact query:

```text
CFTR[gene] AND "conflicting classifications of pathogenicity"[clinsig]
```

The query returned 549 Variation IDs. A search match may reflect a conflicting
condition-level RCV assertion even when the VCV aggregate germline label is not
conflicting. ESummary screening therefore retained only the 407 records for which
`germline_classification.description` exactly equalled
`Conflicting classifications of pathogenicity`.

The analysis then selected the 25 records with the largest counts of supporting
SCV accessions. Versioned VCV accession was the deterministic tie-breaker. This
rule creates an explicitly disclosed genuine subset focused on records with the
largest evidence-reconciliation workload.

## Retrieval procedure

1. `esearch` obtains the total count and complete Variation ID list.
2. `esummary` retrieves JSON records in batches of 100.
3. Summaries are screened for the exact aggregate conflict label.
4. Eligible records are ranked by the number of supporting SCV accessions.
5. `efetch` retrieves one VCV XML archive for each selected Variation ID.
6. Requests use bounded exponential retries, a descriptive tool identifier, a
   user agent, and pacing compatible with NCBI's public service.

The optional `--email` value is sent only as an NCBI contact parameter; it is not
written to an output file.

## Validation rules

### Search and summary validation

- Returned ID count must equal the count reported by ESearch.
- IDs must be unique.
- Every requested ID must have an ESummary response.
- A versioned accession must match `VCV` followed by nine digits and a version.
- CFTR must be present in the gene assignments.
- Aggregate germline classification must exactly match the conflict label.

### Full-record validation

- XML must contain exactly one `VariationArchive`.
- `RecordStatus` must equal `current`; previous or retired records cause failure.
- The record must have a CFTR gene assignment.
- The VCV XML accession version must match ESummary.
- The aggregate germline classification must remain conflicting.
- Classification composition must follow ClinVar's documented
  `Label (count); Label (count)` structure.

### Coordinate handling

The parser selects the `GRCh38` `SequenceLocation` whose `AssemblyStatus` is
`current`, preferring the location marked `forDisplay="true"`. The assembly
accession, one-based VCF position, reference allele, and alternate allele are
reported separately. The canonical SPDI field retains its native interbase
coordinate representation and is not silently substituted for the VCF position.

## Extracted variables

- Variation ID and versioned VCV accession
- record status, variant name, and variant type
- canonical SPDI, protein change, and dbSNP identifiers
- molecular consequence annotations
- current GRCh38 assembly accession and alleles
- aggregate classification, review status, and last-evaluated date
- submission and submitter counts
- grouped classification counts
- condition names, PubMed IDs, and maximum reported allele frequency
- current-record, criteria-provided, polarization, and staleness flags

`not provided` and `not specified` labels are excluded from the informative-trait
count but remain available on the original ClinVar page.

## Classification grouping

ClinVar's aggregate conflict explanation is decomposed into four reporting groups:

| Reporting group | Included labels |
| --- | --- |
| Pathogenic tier | `Pathogenic`, `Likely pathogenic` |
| Uncertain significance | `Uncertain significance` |
| Benign tier | `Benign`, `Likely benign` |
| Other | Any remaining aggregate labels |

`polarized_pathogenic_benign` is true only when at least one pathogenic-tier and
at least one benign-tier classification are represented. This flag describes the
breadth of disagreement; it does not adjudicate the variant.

## Manual-curation priority rules

| Tier | Rule |
| --- | --- |
| Tier 1 | Criteria provided, at least three submitters, and pathogenic-versus-benign polarization |
| Tier 2 | Criteria provided and at least three submitters, without the Tier 1 polarization |
| Tier 3 | More limited review support than Tier 1 or Tier 2 |

Within the source subset, records are ordered by total aggregate submission count
and versioned accession. The tiers are workflow rules designed to organize manual
review; they are not a pathogenicity model, diagnostic score, or clinical
recommendation.

## Verified results

- 549 records matched the live search query.
- 407 had the exact aggregate germline conflict label in ESummary.
- 25 high-submission records were selected and validated as current.
- 25 had current GRCh38 coordinates.
- 25 had criteria-provided review status.
- 7 had pathogenic-versus-benign polarization and entered Tier 1.
- 18 entered Tier 2.
- 0 were more than five years past the aggregate last-evaluated date.
- Median aggregate submission count was 25; maximum was 35.
- Consequences comprised 20 missense, 4 synonymous, and 1 intronic annotation.

## Output generation

The main record table preserves extracted fields and all decision flags. A second,
long-format table records classification composition. The figure shows composition
for the 15 highest-submission records and the complete priority-tier distribution.
Both PNG and SVG formats are written from the same plotting function.

## Automated validation

Tests use a field-preserving reduced fixture from the authentic 2026-09-28 XML
snapshot for `VCV000035835.90` and the committed output tables. The reduced file
contains only the source elements used by the parser and avoids redistributing
unneeded submitter narrative. The tests verify:

- parsing of identifiers, coordinates, review status, and conflict counts;
- pathogenic-versus-benign polarization and priority tier;
- strict rejection of a non-current record status;
- classification explanation parsing and malformed-input rejection;
- source-manifest/output agreement;
- current status and versioned accessions for all committed records; and
- classification-group totals against aggregate submission counts.

## Reproduction commands

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 analyze_clinvar_conflicts.py
python3 -m unittest discover -s tests -v
```

To write results outside the project directory:

```bash
python3 analyze_clinvar_conflicts.py \
  --data-dir /path/to/data \
  --output-dir /path/to/outputs
```

## Reproducibility boundary

ClinVar's live API follows its release cycle. The code reproduces the method, while
the committed CSV files and source manifest preserve the verified 2026-09-28
result. A future live run may select newer VCV versions or a different leading 25.
Exact snapshot comparison should use the accession versions in the manifest.

## Limitations and safeguards

- A high-submission subset emphasizes curation workload and is not representative
  of all CFTR variation.
- Condition-specific assertions may be combined in a VCV aggregate.
- Submission counts are not votes and are not weighted by laboratory performance.
- Reported allele frequencies are descriptive source fields, not a harmonized
  ancestry-aware frequency analysis.
- No assertion is reclassified by this workflow.
- No patient-level information is retrieved or processed.
- Outputs are unsuitable for diagnosis, treatment, or direct clinical reporting.
