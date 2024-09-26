.PHONY: all setup clean

all: setup

setup:
	@echo "Setting up the RNASeq Data Analyser environment..."
	conda env create -f environment.yml
	@mkdir -p data results logs external
	@echo "Setup complete."

clean:
	@echo "Cleaning up..."
	conda remove --name rnaseq_env --all -y
	rm -rf data/* results/* logs/* external/*
	@echo "Clean up complete."
