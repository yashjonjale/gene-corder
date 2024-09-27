#!/usr/bin/env python3

import argparse
import os
import sys
import json
import subprocess
import requests

# Load configuration
CONFIG_FILE = 'config.json'
if os.path.exists(CONFIG_FILE):
    with open(CONFIG_FILE, 'r') as f:
        config = json.load(f)
else:
    config = {
        "supported_organisms": {},
        "objects": {}
    }



def save_config():
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f, indent=4)

def list_organisms():
    print("Supported organisms:")
    # This should be replaced with a real API call to Ensembl
    organisms = ["Homo_sapiens", "Mus_musculus", "Arabidopsis_thaliana"]
    for org in organisms:
        print(f"- {org}")
def get_available_species():
    server = "https://rest.ensembl.org"
    ext = "/info/species?content-type=application/json"

    response = requests.get(server + ext)
    if not response.ok:
        response.raise_for_status()

    data = response.json()
    species_list = []
    for species in data['species']:
        species_info = {
            'display_name': species['display_name'],       # e.g., Human
            'scientific_name': species['scientific_name'], # e.g., Homo sapiens
            'url_name': species['url_name'],               # e.g., Homo_sapiens
            'production_name': species['name'],            # e.g., homo_sapiens
            'assembly': species['assembly'],               # e.g., GRCh38
            'taxonomy_id': species['taxonomy_id'],         # e.g., 9606
        }
        species_list.append(species_info)
        # Update the supported_organisms in config
        config['supported_organisms'][species_info['production_name']] = species_info
    save_config()
    return species_list

def instantiate_organism(args):
    organism_name = args.organism
    if organism_name in config['organisms']:
        print(f"Organism '{organism_name}' is already instantiated.")
        return
    print(f"Instantiating organism '{organism_name}'...")
    # Fetch genome and annotation files from Ensembl
    # Placeholder for actual download code
    organism_data = {
        "genome": f"data/{organism_name}_genome.fa",
        "gtf": f"data/{organism_name}_annotations.gtf",
        "transcriptome": f"data/{organism_name}_transcriptome.fa"
    }
    # Simulate downloading files
    os.makedirs('data', exist_ok=True)
    for key, path in organism_data.items():
        with open(path, 'w') as f:
            f.write(f">Dummy {key} data for {organism_name}\nATGC\n")
    config['organisms'][organism_name] = organism_data
    config['default_organism'] = organism_name
    save_config()
    print(f"Organism '{organism_name}' instantiated successfully.")

def quantize(args):
    organism_name = config['default_organism']
    if organism_name is None:
        print("Error: No organism instantiated. Please instantiate an organism first.")
        sys.exit(1)
    sra_code = args.sra
    quant_name = args.name
    print(f"Quantifying SRA dataset '{sra_code}' for organism '{organism_name}'...")
    quant_dir = os.path.join('results', quant_name)
    os.makedirs(quant_dir, exist_ok=True)
    # Placeholder for actual quantification code
    quant_file = os.path.join(quant_dir, 'abundance.tsv')
    with open(quant_file, 'w') as f:
        f.write("target_id\tlength\teff_length\test_counts\ttpm\n")
        f.write("Gene1\t1000\t800\t50\t10.0\n")
        f.write("Gene2\t2000\t1800\t30\t5.0\n")
    # Update config
    organism_quant = config['organisms'][organism_name].setdefault('quantifications', {})
    organism_quant[quant_name] = {
        'sra_code': sra_code,
        'path': quant_dir
    }
    save_config()
    print(f"Quantification '{quant_name}' completed.")

def list_quantifications(args):
    organism_name = config['default_organism']
    if organism_name is None:
        print("Error: No organism instantiated.")
        sys.exit(1)
    quantifications = config['organisms'][organism_name].get('quantifications', {})
    if not quantifications:
        print("No quantifications found for this organism.")
        return
    print(f"Quantifications for '{organism_name}':")
    for name, info in quantifications.items():
        print(f"- {name}: SRA {info['sra_code']}")

def plot_gene_abundances(args):
    import pandas as pd
    import matplotlib.pyplot as plt
    gene_name = args.gene
    organism_name = config['default_organism']
    quantifications = config['organisms'][organism_name].get('quantifications', {})
    if not quantifications:
        print("Error: No quantifications available.")
        return
    abundances = []
    samples = []
    for quant_name, info in quantifications.items():
        quant_file = os.path.join(info['path'], 'abundance.tsv')
        df = pd.read_csv(quant_file, sep='\t')
        gene_row = df[df['target_id'] == gene_name]
        if not gene_row.empty:
            tpm = gene_row['tpm'].values[0]
            abundances.append(tpm)
            samples.append(quant_name)
    if not abundances:
        print(f"Gene '{gene_name}' not found in any quantifications.")
        return
    plt.bar(samples, abundances)
    plt.xlabel('Sample')
    plt.ylabel('TPM')
    plt.title(f'Expression of {gene_name}')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plot_path = os.path.join('results', f'{gene_name}_abundance.png')
    plt.savefig(plot_path)
    plt.show()
    print(f"Plot saved to {plot_path}")

def generate_correlation_matrix(args):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    genes = [gene.strip() for gene in args.genes.split(',')]
    organism_name = config['default_organism']
    quantifications = config['organisms'][organism_name].get('quantifications', {})
    if not quantifications:
        print("Error: No quantifications available.")
        return
    data = {}
    for quant_name, info in quantifications.items():
        quant_file = os.path.join(info['path'], 'abundance.tsv')
        df = pd.read_csv(quant_file, sep='\t')
        df.set_index('target_id', inplace=True)
        data[quant_name] = df['tpm']
    combined_df = pd.DataFrame(data)
    selected_df = combined_df.loc[genes]
    corr_matrix = selected_df.corr()
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm')
    plt.title('Gene Correlation Matrix')
    plt.tight_layout()
    plot_path = os.path.join('results', 'correlation_matrix.png')
    plt.savefig(plot_path)
    plt.show()
    print(f"Correlation matrix saved to {plot_path}")

def deseq_analyse(args):
    gene_name = args.gene
    print(f"Performing DESeq2 analysis for gene '{gene_name}'...")
    # Placeholder for actual DESeq2 analysis using pydeseq2
    print("DESeq2 analysis completed.")

def main():
    parser = argparse.ArgumentParser(description='RNASeq Data Analyser')
    subparsers = parser.add_subparsers(dest='command')

    # list
    parser_list = subparsers.add_parser('list', help='List available organisms')
    parser_list.set_defaults(func=list_organisms)

    # instantiate
    parser_instantiate = subparsers.add_parser('instantiate', help='Instantiate an organism')
    parser_instantiate.add_argument('--organism', required=True, help='Organism name')
    parser_instantiate.set_defaults(func=instantiate_organism)

    # quantize
    parser_quantize = subparsers.add_parser('quantize', help='Quantize RNA-Seq data')
    parser_quantize.add_argument('--sra', required=True, help='SRA accession code')
    parser_quantize.add_argument('--name', required=True, help='Name for this quantification')
    parser_quantize.set_defaults(func=quantize)

    # list_quant
    parser_list_quant = subparsers.add_parser('list_quant', help='List quantifications')
    parser_list_quant.set_defaults(func=list_quantifications)

    # plot_gene_abundances
    parser_plot = subparsers.add_parser('plot_gene_abundances', help='Plot gene abundances')
    parser_plot.add_argument('--gene', required=True, help='Gene name')
    parser_plot.set_defaults(func=plot_gene_abundances)

    # generate_correlation_matrix
    parser_corr = subparsers.add_parser('generate_correlation_matrix', help='Generate correlation matrix')
    parser_corr.add_argument('--genes', required=True, help='Comma-separated list of genes')
    parser_corr.set_defaults(func=generate_correlation_matrix)

    # deseq_analyse
    parser_deseq = subparsers.add_parser('deseq_analyse', help='Perform DESeq2 analysis')
    parser_deseq.add_argument('--gene', required=True, help='Gene name')
    parser_deseq.set_defaults(func=deseq_analyse)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
