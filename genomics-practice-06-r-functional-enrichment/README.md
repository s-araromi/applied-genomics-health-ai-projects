# Project 06: Functional Enrichment of Paired Lung-Tumor Expression Signals

## Project overview

This project interprets the robust differential-expression signals identified in
Project 05 by testing whether specific Gene Ontology (GO) biological processes and
Reactome pathways occur more often than expected. The analysis treats genes higher
and lower in tumor tissue as separate sets and uses the genes measurable on the
GPL570 microarray as the background.

The implementation uses base R, retrieves authoritative annotation files, validates
their byte counts and checksums, propagates GO annotations through safe ontology
relations, controls the false-discovery rate, and reports redundancy-aware
representative terms.

## Analytical question

Which biological processes and pathways are over-represented among the robust genes
with higher or lower expression in lung tumors relative to paired adjacent-normal
tissue in GSE19804?

## Biological and clinical context

Project 05 identified 1,328 unique, unambiguous genes that met the declared effect
size and false-discovery criteria in both the 60-pair primary analysis and the
53-pair sensitivity analysis. Of these genes, 430 were higher and 898 were lower in
tumor tissue.

Over-representation analysis asks whether an annotation category contains more
input genes than expected under a defined background. It provides a structured
summary of a gene list; it does not directly measure pathway activity or establish
causal disease mechanisms.

## Data source

| Input | Exact record or release | Direct source |
| --- | --- | --- |
| Project 05 robust genes and measured universe | Derived from complete GSE19804 and GPL570 Project 05 results | [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804), [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) |
| Human GO annotations | GOA human GAF generated 20 May 2026; declared GO version `2026-04-27` | [goa_human.gaf.gz](https://current.geneontology.org/annotations/goa_human.gaf.gz) |
| GO ontology | `go-basic.obo`, release `2026-07-26` | [go-basic.obo](https://current.geneontology.org/ontology/go-basic.obo) |
| Reactome mappings | NCBI-to-Reactome all-level mapping, release 97 | [NCBI2Reactome_All_Levels.txt](https://reactome.org/download/97/NCBI2Reactome_All_Levels.txt) |

Two compact Project 05-derived inputs are tracked in `data/`:

- `project05_analysis_universe.csv`: 20,848 unique GPL570 genes eligible for
  statistical testing.
- `project05_robust_gene_set.csv`: 1,328 robust unique genes, including their
  directions, selected probes, effect estimates, and FDR values.

They are deterministic outputs of the genuine GSE19804 analysis, not simulated or
substituted data.

## Reference and data-quality assessment

All source URLs were active and the complete files were downloaded and validated on
21 September 2026. The workflow pins exact byte counts, MD5 values, and SHA-256
values. A source file is rejected if any expected value differs.

The GAF header reports generation on 20 May 2026, an underlying UniProt source-file
date of 30 April 2026, and GO version `2026-04-27`. The separately pinned
`go-basic.obo` release is newer (`2026-07-26`); every retained annotation term was
required to remain active in that ontology release.

| File | Bytes | MD5 |
| --- | ---: | --- |
| `goa_human.gaf.gz` | 15,041,996 | `3d16e91814f541afc3cdb03153ae0f55` |
| `go-basic.obo` | 32,227,785 | `26de643e50af9a1588f9a31062348cc4` |
| `NCBI2Reactome_All_Levels.txt` | 98,024,930 | `8ec3d8afe228764cf3111dc2eac6e838` |

The complete SHA-256 values and validation timestamps are recorded in
`data/source_manifest.csv`. Human GO annotations with a `NOT` qualifier were
excluded. GO terms marked obsolete were not used. Only Homo sapiens Reactome
mappings were retained.

## Data usage

Gene Ontology data are distributed under the
[Creative Commons Attribution 4.0 license](https://geneontology.org/docs/go-citation-policy/),
which requires attribution. Reactome database data and derived files are released
under [CC0](https://reactome.org/license); attribution remains encouraged. NCBI
places no restrictions of its own on molecular database data, while noting that
submitters or countries of origin may assert additional rights; see the
[NCBI data-usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).

GSE19804 contains public, de-identified research data. The analysis makes no attempt
to identify participants and does not use restricted clinical records.

## Methods

1. Validate the Project 05 universe and robust-gene files by row count and checksum.
2. Download and validate the GO annotation, GO ontology, and Reactome mapping files.
3. Retain human GO biological-process annotations and remove negated annotations.
4. Propagate direct GO annotations through `is_a` and `part_of` ancestors in the
   acyclic `go-basic` graph.
5. Retain human Reactome mappings across all levels of its pathway hierarchy.
6. Analyze higher-in-tumor and lower-in-tumor genes separately.
7. Restrict each database-specific background to GPL570-measured genes with an
   annotation in that database.
8. Test terms containing 10–500 background genes and at least three input genes
   using the upper-tail hypergeometric test.
9. Apply Benjamini–Hochberg correction within each database and direction.
10. Repeat GO analysis using experimental evidence codes only.
11. Select up to 12 representative significant terms per result group, skipping a
    candidate when the Jaccard similarity of its overlapping input genes is at
    least 0.5 relative to a previously selected term.
12. Export complete results, summaries, a report, and PNG/SVG figures.

The significance criterion was BH FDR ≤ 0.05. Additional technical detail is
provided in [ANALYSIS_WORKFLOW.md](ANALYSIS_WORKFLOW.md).

## Results

| Database | Direction | Annotated background | Annotated input | Tested terms | Significant terms | Significant with experimental GO support |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| GO Biological Process | Higher in tumor | 14,764 | 371 | 1,223 | 255 | 128 |
| Reactome | Higher in tumor | 10,355 | 267 | 288 | 99 | Not applicable |
| GO Biological Process | Lower in tumor | 14,764 | 772 | 2,391 | 756 | 320 |
| Reactome | Lower in tumor | 10,355 | 552 | 497 | 52 | Not applicable |

Selected representative results were:

| Direction | Database | Representative term | Overlap | Enrichment ratio | BH FDR |
| --- | --- | --- | ---: | ---: | ---: |
| Higher in tumor | GO | Mitotic cell cycle process | 47/451 | 4.15 | 9.84 × 10^-14 |
| Higher in tumor | GO | Extracellular matrix organization | 24/208 | 4.59 | 5.55 × 10^-8 |
| Higher in tumor | Reactome | Collagen degradation | 15/64 | 9.09 | 1.21 × 10^-8 |
| Higher in tumor | Reactome | Cell Cycle Checkpoints | 29/274 | 4.10 | 1.21 × 10^-8 |
| Lower in tumor | GO | Vasculature development | 71/308 | 4.41 | 9.40 × 10^-24 |
| Lower in tumor | GO | Blood circulation | 62/358 | 3.31 | 1.36 × 10^-14 |
| Lower in tumor | Reactome | Developmental Cell Lineages | 27/170 | 2.98 | 1.47 × 10^-4 |
| Lower in tumor | Reactome | GPCR ligand binding | 48/432 | 2.08 | 2.41 × 10^-4 |

The complete 4,399-row result table is retained in the outputs rather than limiting
the record to the representative terms.

## Interpretation

Genes higher in tumor were enriched for cell-cycle control, mitotic checkpoints,
collagen remodeling, and extracellular-matrix processes. Genes lower in tumor were
enriched for vascular development, circulation, endothelial and migratory
processes, GPCR signaling, platelet biology, immune signaling, and surfactant
metabolism.

This pattern is compatible with increased proliferation and matrix remodeling in
tumor tissue alongside reduced representation of normal lung vascular and
specialized tissue programs. It may also reflect differences in cell-type
composition between bulk tumor and adjacent-normal samples. It does not establish
that the reported pathways are activated, inhibited, or causal.

## Output files

| File | Description |
| --- | --- |
| `outputs/functional_enrichment_results.csv` | Complete GO and Reactome statistical results |
| `outputs/representative_enriched_terms.csv` | Redundancy-aware representative significant terms |
| `outputs/enrichment_summary.csv` | Background, input, tested-term, and significant-term counts |
| `outputs/functional_enrichment_report.md` | Reproducible results report |
| `outputs/functional_enrichment_overview.png` | Four-panel raster visualization |
| `outputs/functional_enrichment_overview.svg` | Four-panel vector visualization |
| `data/source_manifest.csv` | Source versions, URLs, sizes, checksums, and validation time |

## Reproducibility

### Requirements

- R 4.3 or newer.
- Internet access for the first run.
- `shasum` or `sha256sum` for SHA-256 validation.
- No external R packages.

### Run the analysis

```bash
Rscript analyze_functional_enrichment.R
```

Validated source files are reused on later runs. To force a fresh download:

```bash
Rscript analyze_functional_enrichment.R --refresh-source
```

### Run automated tests

```bash
Rscript tests/test_functional_enrichment.R
```

The test fixtures are documented excerpts from the same GO and Reactome releases
used by the full analysis.

## Technical competencies demonstrated

- Integration of differential-expression results with functional knowledgebases.
- Construction of a platform-aware statistical background.
- Parsing of GAF, OBO, and tab-delimited pathway-mapping formats.
- Traversal of an ontology hierarchy and annotation propagation.
- Direction-stratified hypergeometric testing and BH FDR correction.
- Sensitivity analysis restricted to experimental GO evidence codes.
- Redundancy control using gene-overlap Jaccard similarity.
- Cryptographic source validation and provenance recording.
- Reusable base-R functions, automated testing, and reproducible reporting.
- Cautious interpretation of bulk-tissue expression and pathway annotations.

## Limitations

- The analysis is derived from one microarray cohort of never-smoking women with
  lung adenocarcinoma and is not an independent validation study.
- Over-representation discards gene-level effect-size magnitude after assigning
  genes to a thresholded set.
- GO and Reactome categories overlap and are not statistically independent; the
  representative-term procedure reduces presentation redundancy but does not alter
  the complete hypothesis set.
- Annotation coverage and research attention vary across genes and processes.
- Experimental GO evidence can originate from biological systems other than the
  GSE19804 tissue context.
- Bulk tissue mixes tumor, stromal, immune, vascular, and other cell populations.
- The results do not support individual diagnosis, prognosis, or treatment choice.

## References

- [GSE19804 GEO Series record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804)
- [GPL570 platform record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
- [Gene Ontology download and citation policy](https://geneontology.org/docs/go-citation-policy/)
- [Gene Ontology annotation downloads](https://geneontology.org/docs/download-go-annotations/)
- [Gene Ontology `go-basic` documentation](https://geneontology.org/docs/download-ontology/)
- [Reactome downloads](https://reactome.org/download-data)
- [Reactome license](https://reactome.org/license)
- [NCBI data-usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [Lu et al., 2010, PubMed 20802022](https://pubmed.ncbi.nlm.nih.gov/20802022/)
