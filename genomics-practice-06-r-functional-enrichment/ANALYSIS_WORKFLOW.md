# Analysis Workflow: Direction-Specific GO and Reactome Enrichment

## 1. Analytical objective

Identify GO biological processes and Reactome pathways that are represented more
often than expected among the robust Project 05 genes with higher or lower
expression in GSE19804 lung tumors.

## 2. Data flow

```text
GSE19804 expression matrix + GPL570 annotation
                    |
                    v
        Project 05 robust gene results
        430 higher + 898 lower in tumor
                    |
          +---------+---------+
          |                   |
          v                   v
    GOA + go-basic       Reactome mapping
          |                   |
          +---------+---------+
                    |
                    v
    Direction-specific over-representation
                    |
                    v
 Complete results + evidence sensitivity +
 representative terms + report + visualization
```

## 3. Project 05-derived inputs

The input set contains one selected probe for each robust, uniquely and
unambiguously annotated gene from Project 05. A gene had to meet BH FDR ≤ 0.05 and
absolute mean paired log2 fold change ≥ 1 in both the 60-pair primary analysis and
the 53-pair sensitivity analysis, with the same direction in both analyses.

| Input | Rows | SHA-256 |
| --- | ---: | --- |
| `project05_analysis_universe.csv` | 20,848 | `75eb4546bb2d11a89dfbd2015b2764dea85109a14cc7f51a411c1a5e2176fa11` |
| `project05_robust_gene_set.csv` | 1,328 | `2ff57aa27a8c74516b055cfd725ceb04bd1c38f2dddc52cca9a3ca2dfb439882` |

The analysis is self-contained: these validated derived files are tracked with the
project, so Project 06 does not depend on the local path or runtime outputs of
Project 05.

## 4. Authoritative annotation sources

The analysis uses the GOA human GAF generated on 20 May 2026 (declared GO version
`2026-04-27`), `go-basic.obo` release `2026-07-26`, and the NCBI-to-Reactome
all-level mapping from Reactome release 97, released on 23 June 2026. Every retained
GAF term was required to remain active in the newer pinned OBO release.

Each file is downloaded to `data/` only when absent or when
`--refresh-source` is supplied. A temporary `.download` file is validated before
being renamed into place. Validation includes byte count, MD5, and SHA-256.

The large annotation files are excluded from version control because they total
approximately 145 MB and can be reconstructed from their documented source URLs.

## 5. Gene Ontology processing

The OBO parser retains non-obsolete terms, names, namespaces, and parent relations.
Only biological-process terms are eligible. The GAF parser retains human
(`taxon:9606`) biological-process annotations and discards entries carrying the
`NOT` qualifier.

Direct annotations are propagated to all reachable ancestors through `is_a` and
`part_of` relations. The GO Consortium describes `go-basic` as acyclic and safe for
upward propagation. A cycle check remains in the implementation as an additional
validation guard.

Duplicate gene–term–evidence combinations are removed after propagation.

## 6. Reactome processing

The all-level NCBI mapping is filtered to `Homo sapiens`. Duplicate identifier,
pathway, name, and evidence combinations are removed. NCBI Gene IDs connect the
Project 05 platform annotation to stable Reactome pathway identifiers.

Using the all-level file allows the analysis to test both specific pathways and
their higher-level Reactome categories.

## 7. Statistical background

The background is not the complete human genome. It is the intersection between:

1. the 20,848 unique genes measured and unambiguously annotated on GPL570; and
2. genes represented in the database being tested.

This produces 14,764 GO-annotated background genes and 10,355 Reactome-annotated
background genes. Database-specific input counts are consequently smaller than the
full robust set.

This design reduces selection bias that would result from treating genes that were
not measurable on the array as eligible discoveries.

## 8. Over-representation test

For each direction and database, the one-sided hypergeometric probability is:

```text
P(X >= k) = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
```

where:

- `N` is the annotated background size;
- `K` is the number of background genes assigned to the term;
- `n` is the number of annotated input genes in the direction;
- `k` is the observed input-gene overlap.

Terms were tested when their annotated background size was 10–500 and their input
overlap was at least three. P-values were adjusted with the Benjamini–Hochberg
procedure separately for each database and direction. FDR ≤ 0.05 was significant.

The enrichment ratio is:

```text
(k / n) / (K / N)
```

## 9. Evidence-code sensitivity analysis

The GO analysis was repeated using only the following experimental evidence codes:

```text
EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP
```

A primary significant GO term was labeled experimentally supported when the
corresponding experimental-only result also had FDR ≤ 0.05. This is a robustness
assessment of annotation evidence, not experimental validation within GSE19804.

## 10. Representative-term selection

The complete statistical table retains every tested term. For readable reporting,
significant terms are processed in increasing FDR order. A term is skipped from the
representative set when the Jaccard similarity of its overlapping input genes with
an already selected term is at least 0.5:

```text
J(A, B) = |A intersection B| / |A union B|
```

Up to 12 terms are retained for each database and expression direction. This
procedure changes only the summary and visualization, not significance values.

## 11. Verified analytical results

The full analysis produced 4,399 term-level rows and 48 representative terms.

| Database | Direction | Tested | Significant | Experimental GO support |
| --- | --- | ---: | ---: | ---: |
| GO Biological Process | Higher in tumor | 1,223 | 255 | 128 |
| Reactome | Higher in tumor | 288 | 99 | Not applicable |
| GO Biological Process | Lower in tumor | 2,391 | 756 | 320 |
| Reactome | Lower in tumor | 497 | 52 | Not applicable |

The strongest higher-in-tumor themes included mitotic cell-cycle processes,
checkpoint control, collagen degradation, and matrix organization. The strongest
lower-in-tumor themes included vasculature development, circulation, endothelial
programs, GPCR signaling, immune-related signaling, and surfactant metabolism.

## 12. Output generation

The analysis writes complete and representative result tables, a compact summary,
a Markdown report, a four-panel PNG, and a matching SVG. Point position in the plot
is `-log10(BH FDR)`; point size reflects overlapping-gene count. Orange denotes
higher-in-tumor sets and blue denotes lower-in-tumor sets.

## 13. Automated validation

`tests/test_functional_enrichment.R` contains 11 checks covering:

- pinned source metadata and derived-input checksums;
- the 20,848-gene universe and 1,328-gene robust set;
- OBO hierarchy parsing;
- GAF taxon, aspect, and `NOT` handling;
- `is_a` and `part_of` propagation;
- Reactome human mappings;
- annotated-background hypergeometric calculations;
- GO experimental-evidence sensitivity joins;
- Jaccard similarity and redundancy filtering; and
- final row counts and result summaries.

Fixtures are compact, documented excerpts of the exact authoritative GO and
Reactome source files, not synthetic biological records.

## 14. Reproduction

From the Project 06 directory:

```bash
Rscript analyze_functional_enrichment.R
Rscript tests/test_functional_enrichment.R
```

To redownload and revalidate the source files:

```bash
Rscript analyze_functional_enrichment.R --refresh-source
```

## 15. Interpretation boundaries

Enrichment is conditional on the input thresholds, measured background, annotation
coverage, ontology release, pathway release, and term-size filters. Related terms
share genes, while bulk-tissue expression combines signals from multiple cell
types. These results prioritize biological themes for further study; they do not
constitute mechanistic proof or clinical evidence.

## References

- [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804)
- [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570)
- [GO annotation downloads](https://geneontology.org/docs/download-go-annotations/)
- [GO ontology downloads](https://geneontology.org/docs/download-ontology/)
- [GO citation and license policy](https://geneontology.org/docs/go-citation-policy/)
- [Reactome downloads](https://reactome.org/download-data)
- [Reactome license](https://reactome.org/license)
