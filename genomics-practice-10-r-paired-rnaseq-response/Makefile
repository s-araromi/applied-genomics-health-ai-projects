.PHONY: install test analysis retrieve verify
install:
	Rscript --vanilla install_dependencies.R
test:
	Rscript --vanilla tests/test_rnaseq.R
analysis:
	Rscript --vanilla analyze_rnaseq.R
retrieve:
	Rscript --vanilla retrieve_data.R
verify:
	Rscript --vanilla verify_outputs.R
