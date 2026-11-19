# Data usage and provenance

## Public records

This project uses the complete processed expression matrix and deposited
clinical table for NCBI GEO accession [GSE2034](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034),
plus the [GPL96](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96)
platform annotation. The 286 samples are public GEO accessions. No attempt was
made to identify participants or connect these records to outside information.

The exact download URLs, retrieval date, byte counts, and SHA-256 identities are
recorded in `data/source_manifest.csv`. The full raw downloads are reproducible
with `retrieve_data.py` but are stored in the ignored `data/cache/` directory.
The repository retains a genuine, derived analytical matrix and its hashes so
the verified analysis can run without a network request.

## Reuse terms

The [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)
states that NCBI places no restrictions on GEO data use or distribution while
noting that submitters may retain rights in portions of deposited material.
Accordingly, this project cites the accession and original publication and does
not assert a Creative Commons license that GEO has not granted.

## Outcome handling

GSE2034's public summary reports 180 relapse-free patients and 106 patients who
developed distant metastasis. Its deposited patient table contains 286 rows but
encodes 107 values of `relapse = 1`. This discrepancy is preserved in
`data/reference_assessment.csv`. The analysis uses the row-level `relapse` and
follow-up fields and calls the outcome relapse, not distant metastasis.

For the prespecified five-year classification endpoint, relapse at or before 60
months is positive. A participant followed for at least 60 months without an
earlier relapse is negative at that landmark, including participants who
relapsed later. Ten relapse-free participants censored before 60 months are
excluded because their five-year status is unknown.
