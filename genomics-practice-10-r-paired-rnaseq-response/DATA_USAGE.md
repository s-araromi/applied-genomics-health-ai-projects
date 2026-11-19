# Source attribution and data usage

The analytical input tables are faithful exports of the **airway 1.32.0**
Bioconductor data object, authored and maintained by Michael Love. The original
package DESCRIPTION, including its **LGPL** declaration, is preserved in
`data/airway_DESCRIPTION`. The package's declaration does not specify a
version in that field; this project does not relabel the upstream data under
a different license. The exact source archive and reconstruction code are
linked and pinned in `data/archive_pin.json` and `retrieve_data.R`.

- Package: https://bioconductor.org/packages/airway/
- Versioned source: https://bioconductor.posit.co/packages/3.23/data/experiment/src/contrib/airway_1.32.0.tar.gz
- License information: https://www.gnu.org/licenses/lgpl.html
- Original experiment: Himes BE et al. PLoS ONE 9:e99625 (2014),
  https://doi.org/10.1371/journal.pone.0099625.
- Public record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
- Original reads: https://www.ncbi.nlm.nih.gov/sra?term=SRP033351

NCBI's GEO disclaimer states that NCBI imposes no restrictions on use or
distribution of GEO data, while submitter-held patent, copyright, or other
rights may still apply. NCBI does not independently certify submitted data's
quality or biological interpretation. See:
https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html.

This project redistributes the genuine count and annotation exports with
upstream attribution, original package metadata, exact source links, and
reconstruction instructions. It uses only public donor cell-line labels and
archive accessions. It neither attempts re-identification nor links those
labels to additional personal records. Public availability is not a claim of
unrestricted commercial rights, fresh human-subject consent, or clinical
fitness. The paper's publication license and the package's license describe
different materials; neither is silently applied to everything in GEO.

The live series status was assessed on 19 October 2026. No controlled-access
or confidential patient records were requested, copied, or generated.
