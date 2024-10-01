#!/bin/bash
set -eox pipefail
wget -c "https://raw.githubusercontent.com/yashjonjale/gene-corder/refs/heads/main/env.yml"
eval "$(conda shell.bash hook)"
conda deactivate
conda env create -f env.yml
echo "GeneCorder environment created successfully"
conda activate genecorder
echo "Installing GeneCorder"
pip install git+https://github.com/yashjonjale/gene-corder
echo "To work with the tool activate the environment, run 'conda activate genecorder'"
