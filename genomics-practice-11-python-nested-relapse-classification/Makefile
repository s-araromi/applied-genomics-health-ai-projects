.PHONY: all analyze test verify retrieve

all: analyze test verify

analyze:
	python3 analyze_breast_relapse.py

test:
	python3 -m unittest discover -s tests -v

verify:
	python3 verify_outputs.py

retrieve:
	python3 retrieve_data.py
