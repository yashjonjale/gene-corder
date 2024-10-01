wget -c "https://raw.githubusercontent.com/yashjonjale/gene-corder/refs/heads/main/environment.yml"
eval "$(conda shell.bash hook)"
conda deactivate
conda env create -f environment.yml
echo "GeneCorder environment created successfully"
conda activate gc
echo "Installing GeneCorder"
pip install git+https://github.com/yashjonjale/gene-corder
echo "To work with the tool activate the environment, run 'conda activate gc'"
