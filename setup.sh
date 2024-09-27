#!/bin/bash

eval "$(conda shell.bash hook)"

conda deactivate
conda env create -f environment.yml
conda activate $1


