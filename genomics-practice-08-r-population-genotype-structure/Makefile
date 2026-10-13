.PHONY: all verify test analyze refresh

all: verify test analyze

verify:
	python3 retrieve_source_subset.py

test:
	Rscript --vanilla tests/test_genotype_functions.R

analyze:
	Rscript --vanilla analyze_population_structure.R

# Network retrieval is explicit; the normal analysis uses the pinned local data.
refresh:
	python3 retrieve_source_subset.py --refresh
