.PHONY: all analyze test verify

all: analyze test verify

analyze:
	python3 analyze_cross_cohort_replication.py

test:
	python3 -m unittest discover -s tests -v

verify:
	python3 analyze_cross_cohort_replication.py --no-download
