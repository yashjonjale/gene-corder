.PHONY: clean all setup

ENV_NAME = "gc"
all: setup
	@echo "Building..."

setup:
	@echo "Setting up..."
	@bash setup.sh $(ENV_NAME)
	python setup.py 
	@echo "Done. Activate the environment with: conda activate $(ENV_NAME)"
