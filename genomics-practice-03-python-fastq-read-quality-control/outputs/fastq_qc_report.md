# Paired-End FASTQ Quality-Control Report

## Dataset

- ENA run: `ERR194147`
- Sample: NA12878
- Subset: first 5,000 read pairs
- Read length: 101–101 bases (R1); 101–101 bases (R2)

## Results

| Mate | Reads | Bases | Mean Phred | Q30 bases | GC content | N bases | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| R1 | 5,000 | 505,000 | 35.93 | 92.80% | 49.87% | 12 | PASS |
| R2 | 5,000 | 505,000 | 35.53 | 91.53% | 49.74% | 1,292 | REVIEW |
| Combined | 10,000 | 1,010,000 | 35.73 | 92.16% | 49.81% | 1,304 | REVIEW |

## Interpretation

Q30 represents an estimated base-call error probability of 0.1%. A mate is marked
`REVIEW` when fewer than 90% of bases reach Q30 or more than 0.1% are ambiguous
`N` calls. A review flag is descriptive and does not by itself make the full run
unsuitable for downstream analysis.
