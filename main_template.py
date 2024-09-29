#!/usr/bin/env python3

import argparse
import os
import sys
import json
import subprocess
import requests
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pytximport as txi


def save_config(config_json):
    ##save the config file
    name = 'config.json'
    #delete the existing file config.json
    if os.path.exists(name):
        os.remove(name)
    with open(name, 'w') as f:
        json.dump(config_json, f, indent=4)
    






def instantiate_organism(args):
    print(f"Instantiating Organism ....")
    
    organism_name = args.organism
    name = args.name
    desc = args.desc
    transcript_path = args.transcriptome_path
    genome_path = args.genome_path
    gtf_path = args.gtf_path

    print(f"Organism: {organism_name}")
    print(f"Name: {name}")
    print(f"Instantiating object {name} for organism {organism_name}...")
    #debug output
    print(f"Transcriptome Path: {transcript_path}")
    print(f"Genome Path: {genome_path}")
    print(f"GTF Path: {gtf_path}")

    # Check if paths are global, if not, make them global
    if transcript_path is not None:
        if not os.path.isabs(transcript_path):
            transcript_path = os.path.abspath(transcript_path)
    if not os.path.isabs(genome_path):
        genome_path = os.path.abspath(genome_path)
   
    if not os.path.isabs(gtf_path):
        gtf_path = os.path.abspath(gtf_path)

    # Create object for config
    obj = {
        "name": name,
        "desc": desc,
        "organism": organism_name,
        "paths": {
            "transcriptome": transcript_path,
            "genome": genome_path,
            "gtf": gtf_path
        },
        "alt_paths": {
            "transcriptome": None,
            "genome": None,
            "gtf": None
        },
        "quantifications": {
            # "quant1": {
            #     "sra_codes": [{
            #         "sra_code": "SRR12345",
            #         "path": "path/to/quant1"
                    #   "path_to_fastq": ,
            #     }, {
            #         "sra_code": "SRR67890",
            #         "name": "Sample2"
            #     }],
            #     "type": "paired",
            # }
        },
        "index_paths": {
            "kallisto_index": None,
            "t2g": None
        }
    }

    # Load the current config and append the new object
    config = None
    with open('config.json', 'r') as f:
        config = json.load(f)
        # config['objects'].append(obj)
    config['objects'][name] = obj
    # Directory to store outputs
    output_dir = f"./data/{name}"
    os.makedirs(output_dir, exist_ok=True)

    ## gffread -w ehuxleyi_transcripts.fa -g Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf
    ## run the above using subprocess for the obj paths
    if transcript_path is None or not os.path.exists(transcript_path):
        print(f"Converting GTF to FASTA...")
        gffread_cmd = [
            "gffread", "-w", f"{output_dir}/{name}_transcripts.fa", "-g", genome_path, gtf_path
        ]
        try:
            subprocess.run(gffread_cmd, check=True)
            print(f"FASTA file generated at {output_dir}/{name}_transcripts.fa")
            obj["paths"]["transcriptome"] = f"{output_dir}/{name}_transcripts.fa"
            transcript_path = f"{output_dir}/{name}_transcripts.fa"
            #make global path
            transcript_path = os.path.abspath(transcript_path)
            print(f"Transcript FastA path: {transcript_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running gffread: {e}")
            return
    

    # Step 1: Index genome using kallisto index
    kallisto_index_path = os.path.join(output_dir, f"index.idx")
    #make global path
    kallisto_index_path = os.path.abspath(kallisto_index_path)
    print(f"Indexing genome with kallisto...")
    kallisto_cmd = [
        "kallisto", "index",
        "-i", kallisto_index_path,
        transcript_path
    ]
    
    try:
        subprocess.run(kallisto_cmd, check=True)
        print(f"Kallisto index created at {kallisto_index_path}")
        obj["index_paths"]["kallisto_index"] = kallisto_index_path
    except subprocess.CalledProcessError as e:
        print(f"Error while running kallisto index: {e}")
        return

    # Step 2: Generate t2g.tsv using kb ref
    t2g_path = os.path.join(output_dir, f"{name}_t2g.tsv")
    #make global path
    t2g_path = os.path.abspath(t2g_path)
    print(f"Generating t2g.tsv using kb ref...")
    kb_ref_cmd = [
        "kb", "ref", "-i", "ind.idx",
        "-g", t2g_path,
        "-f1", transcript_path,
        genome_path, gtf_path
    ]
    
    try:
        subprocess.run(kb_ref_cmd, check=True)
        print(f"t2g.tsv generated at {t2g_path}")
        obj["index_paths"]["t2g"] = t2g_path
        #delete ind.idx
        os.remove("ind.idx")
    except subprocess.CalledProcessError as e:
        print(f"Error while running kb ref: {e}")
        return

    # Step 3: Save the updated config file with the new paths
    #debug print obj
    print("Updated object:")
    print(obj)
    config['objects'][name] = obj
    save_config(config)
    print(f"Updated config file with new paths for {name}.")
    
    print(f"Genome indexing and t2g generation complete. Ready for downstream quantification!")


def quantize(args):
    sra_list = args.sra.split(',') 
    obj_name = args.obj
    name = args.name
    paired = args.paired
    print("Quantizing RNA-Seq data...")
    
    ##debug output
    print(f"Name: {name}")
    print(f"SRA List: {sra_list}")
    print(f"Object Name: {obj_name}")

    config = None
    with open('config.json', 'r') as f:
        config = json.load(f)

    obj = config['objects'][obj_name]

    ## Object is fetched. 
    ## You have to download the SRA files into the data folder and save the path in the obj quantification path entries

    for sra in sra_list:
        #create path for the sra file, store it in ./data/<obj_name>/<name>/<sra_code>.sra
        sra_path = f"./data/{obj_name}/{name}/"
        os.makedirs(os.path.dirname(sra_path), exist_ok=True)
        
        download = 0  #Downloading required ?
        
        #create folder if it does not exist
        if not os.path.exists(os.path.dirname(sra_path + f"{sra}/")):
            download = 1 #download required 
            
        #download the sra file using subprocess and prefetch
        if download == 1:
            #Downloading
            print("Downlading ", sra)

            prefetch_cmd = [
                "prefetch",
                sra,
                "-O",
                sra_path
            ]

            try:
                subprocess.run(prefetch_cmd, check=True)
                print(f"SRA file {sra} downloaded to {sra_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error while downloading SRA file {sra}: {e}")
                return
    
     ## once downloaded in the corresponding folder, you have to do fastqdump on the corresponding sra files using subprocess

    for sra in sra_list:
        sra_path = f"./data/{obj_name}/{name}/{sra}/{sra}.sra"
        fastqdump_path = f"./data/{obj_name}/{name}/{sra}/"
        
        if paired:
            fastqdump_cmd = [
                "fastq-dump",
                "--split-files",
                "-O",
                fastqdump_path,
                sra_path
            ]
        else:
            fastqdump_cmd = [
                "fastq-dump",
                "-O",
                fastqdump_path,
                sra_path
            ]
        
        try:
            subprocess.run(fastqdump_cmd, check=True)
            print(f"Fastq files generated for {sra}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running fastq-dump for {sra}: {e}")
            return
        
    ## Now, quantfying the fastq files using kallisto, different commands for paired and single end reads

    ## you have to save the path to the quantification in the obj quantification path entries

    ## run kallisto quantification for each fastq file
    abundances_tsv_paths = []
    for sra in sra_list:
        fastqdump_path = f"./data/{obj_name}/{name}/{sra}/"
        kallisto_index_path = obj["index_paths"]["kallisto_index"]
        output_dir = f"./data/{obj_name}/{name}/{sra}_quant/"
        os.makedirs(output_dir, exist_ok=True)
        if paired:
            kallisto_cmd = [
                "kallisto", "quant",
                "-i", kallisto_index_path,
                "-o", output_dir,
                f"{fastqdump_path}/{sra}_1.fastq",
                f"{fastqdump_path}/{sra}_2.fastq"
            ]
        else:
            kallisto_cmd = [
                "kallisto", "quant",
                "-i", kallisto_index_path,
                "-o", output_dir,
                f"{fastqdump_path}/{sra}.fastq"
            ]
        
        try:
            subprocess.run(kallisto_cmd, check=True)
            print(f"Quantification completed for {sra}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running kallisto quant for {sra}: {e}")
            return

        ## save the path to the quantification in the obj quantification path entries
        ## make the new entry in the quantifications dictionary accordingly

        obj["quantifications"][name]["paths"].append(output_dir)
        if paired:
            obj["quantifications"][name]["types"].append("paired")
        else:
            obj["quantifications"][name]["types"].append("single")

        abundances_tsv_path = os.path.join(output_dir, "abundance.tsv")
        #make global path
        abundances_tsv_path = os.path.abspath(abundances_tsv_path)
        abundances_tsv_paths.append(abundances_tsv_path)
    
        

def run_liftoff(ref_genome, ref_annotations, target_genome, output_annotations, threads=1, min_coverage=0.95, copies=False):
    # Base liftoff command
    liftoff_cmd = [
        'liftoff',
        '-g', ref_annotations,            # Path to the reference annotations (GFF/GTF)
        '-o', output_annotations,         # Path to the output annotations file (GFF/GTF)
        '-p', str(threads),               # Number of threads
        '-e', str(min_coverage),          # Minimum coverage percentage of gene sequence in target genome
        target_genome,                    # Target genome file (FASTA)
        ref_genome                        # Reference genome file (FASTA)
    ]
    
    # Option to include the --copies flag if we want to map multiple copies of the same gene
    if copies:
        liftoff_cmd.append('--copies')

    # Run the liftoff command
    try:
        result = subprocess.run(liftoff_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"Liftoff completed successfully!\nOutput:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Liftoff:\n{e.stderr}")

# Example usage
ref_genome = "ref_genome.fasta"
ref_annotations = "ref_genome.gff"
target_genome = "target_genome.fasta"
output_annotations = "target_genome_annotations.gff"

# Call the function to run Liftoff
run_liftoff(ref_genome, ref_annotations, target_genome, output_annotations, threads=8, min_coverage=0.95, copies=True)

    

def intersect(args):
    #capture the arguments
    alt_genome_path = args.alt_genome_path
    alt_gtf_path = args.alt_gtf_path
    obj_name = args.obj
    #loaf the config file
    config = None
    with open('config.json', 'r') as f:
        config = json.load(f)
    
    #fetch the object
    obj = config['objects'][obj_name]
    #make paths global
    if not os.path.isabs(alt_genome_path):
        alt_genome_path = os.path.abspath(alt_genome_path)
    if not os.path.isabs(alt_gtf_path):
        alt_gtf_path = os.path.abspath(alt_gtf_path)
    
    
    









    # organism_name = config['default_organism']
    # if organism_name is None:
    #     print("Error: No organism instantiated. Please instantiate an organism first.")
    #     sys.exit(1)
    # sra_code = args.sra
    # quant_name = args.name
    # print(f"Quantifying SRA dataset '{sra_code}' for organism '{organism_name}'...")
    # quant_dir = os.path.join('results', quant_name)
    # os.makedirs(quant_dir, exist_ok=True)
    # # Placeholder for actual quantification code
    # quant_file = os.path.join(quant_dir, 'abundance.tsv')
    # with open(quant_file, 'w') as f:
    #     f.write("target_id\tlength\teff_length\test_counts\ttpm\n")
    #     f.write("Gene1\t1000\t800\t50\t10.0\n")
    #     f.write("Gene2\t2000\t1800\t30\t5.0\n")
    # # Update config
    # organism_quant = config['organisms'][organism_name].setdefault('quantifications', {})
    # organism_quant[quant_name] = {
    #     'sra_code': sra_code,
    #     'path': quant_dir
    # }
    # save_config()
    # print(f"Quantification '{quant_name}' completed.")

def list_quantifications(args):
    print("Listing quantifications...")
    # organism_name = config['default_organism']
    # if organism_name is None:
    #     print("Error: No organism instantiated.")
    #     sys.exit(1)
    # quantifications = config['organisms'][organism_name].get('quantifications', {})
    # if not quantifications:
    #     print("No quantifications found for this organism.")
    #     return
    # print(f"Quantifications for '{organism_name}':")
    # for name, info in quantifications.items():
    #     print(f"- {name}: SRA {info['sra_code']}")

def plot_gene_abundances(args):
    print("Plotting gene abundances...")
    # import pandas as pd
    # import matplotlib.pyplot as plt
    # gene_name = args.gene
    # organism_name = config['default_organism']
    # quantifications = config['organisms'][organism_name].get('quantifications', {})
    # if not quantifications:
    #     print("Error: No quantifications available.")
    #     return
    # abundances = []
    # samples = []
    # for quant_name, info in quantifications.items():
    #     quant_file = os.path.join(info['path'], 'abundance.tsv')
    #     df = pd.read_csv(quant_file, sep='\t')
    #     gene_row = df[df['target_id'] == gene_name]
    #     if not gene_row.empty:
    #         tpm = gene_row['tpm'].values[0]
    #         abundances.append(tpm)
    #         samples.append(quant_name)
    # if not abundances:
    #     print(f"Gene '{gene_name}' not found in any quantifications.")
    #     return
    # plt.bar(samples, abundances)
    # plt.xlabel('Sample')
    # plt.ylabel('TPM')
    # plt.title(f'Expression of {gene_name}')
    # plt.xticks(rotation=45)
    # plt.tight_layout()
    # plot_path = os.path.join('results', f'{gene_name}_abundance.png')
    # plt.savefig(plot_path)
    # plt.show()
    # print(f"Plot saved to {plot_path}")

def generate_correlation_matrix(args):\
    print("Generating correlation matrix...")
    
    # genes = [gene.strip() for gene in args.genes.split(',')]
    # organism_name = config['default_organism']
    # quantifications = config['organisms'][organism_name].get('quantifications', {})
    # if not quantifications:
    #     print("Error: No quantifications available.")
    #     return
    # data = {}
    # for quant_name, info in quantifications.items():
    #     quant_file = os.path.join(info['path'], 'abundance.tsv')
    #     df = pd.read_csv(quant_file, sep='\t')
    #     df.set_index('target_id', inplace=True)
    #     data[quant_name] = df['tpm']
    # combined_df = pd.DataFrame(data)
    # selected_df = combined_df.loc[genes]
    # corr_matrix = selected_df.corr()
    # plt.figure(figsize=(8, 6))
    # sns.heatmap(corr_matrix, annot=True, cmap='coolwarm')
    # plt.title('Gene Correlation Matrix')
    # plt.tight_layout()
    # plot_path = os.path.join('results', 'correlation_matrix.png')
    # plt.savefig(plot_path)
    # plt.show()
    # print(f"Correlation matrix saved to {plot_path}")

def deseq_analyse(args):
    print("Performing DESeq2 analysis...")
    # gene_name = args.gene
    # print(f"Performing DESeq2 analysis for gene '{gene_name}'...")
    # # Placeholder for actual DESeq2 analysis using pydeseq2
    # print("DESeq2 analysis completed.")

def main():

    #check if config file exists
    if not os.path.exists('config.json'):
        with open('config.json', 'w') as f:
            json.dump({"objects": {}}, f, indent=4)
    

    parser = argparse.ArgumentParser(description='RNASeq Data Analyser')
    subparsers = parser.add_subparsers(dest='command')

    # list - removed
    # parser_list = subparsers.add_parser('list', help='List available organisms')
    # parser_list.set_defaults(func=list_organisms)

    # instantiate
    parser_instantiate = subparsers.add_parser('instantiate', help='Instantiate an object for an organism')
    parser_instantiate.add_argument('--organism', required=True, help='Organism name')
    parser_instantiate.add_argument('--name', required=True, help='Name for this quantification')
    parser_instantiate.add_argument('--desc', required=False, help='Description for this quantification')
    parser_instantiate.add_argument('--transcriptome_path', required=True, help='path to transcriptome file in fasta format')
    parser_instantiate.add_argument('--genome_path', required=True, help='path to genome file in fasta format')
    parser_instantiate.add_argument('--gtf_path', required=True, help='path to the annotation file in gtf/gff format')
    parser_instantiate.set_defaults(func=instantiate_organism)

    # quantize
    parser_quantize = subparsers.add_parser('quantize', help='Quantize RNA-Seq data')
    parser_quantize.add_argument('--sra', required=False, help='SRA accession code')
    parser_quantize.add_argument('--name', required=False, help='Name for this quantification')
    parser_quantize.add_argument('--obj', required=False, help='Object name')
    parser_quantize.add_argument('--paired', required=False, action='store_true', help='Paired-end reads')

    parser_quantize.set_defaults(func=quantize)

    # # list_quant
    # parser_list_quant = subparsers.add_parser('list_quant', help='List quantifications')
    # parser_list_quant.set_defaults(func=list_quantifications)



    # # intersect
    parser_intersect = subparsers.add_parser('intersect', help='Intersect gene sets')
    parser_intersect.add_argument('--alt_genome_path', required=True, help='Alternative genome path')
    parser_intersect.add_argument('--alt_gtf_path', required=True, help='Alternative GTF path')
    parser_intersect.add_argument('--obj', required=True, help='Object name')
    parser_intersect.set_defaults(func=intersect)

    # # plot_gene_abundances
    # parser_plot = subparsers.add_parser('plot_gene_abundances', help='Plot gene abundances')
    # parser_plot.add_argument('--gene', required=False, help='Gene name')
    # parser_plot.set_defaults(func=plot_gene_abundances)

    # # generate_correlation_matrix
    # parser_corr = subparsers.add_parser('generate_correlation_matrix', help='Generate correlation matrix')
    # parser_corr.add_argument('--genes', required=False, help='Comma-separated list of genes')
    # parser_corr.set_defaults(func=generate_correlation_matrix)

    # # deseq_analyse
    # parser_deseq = subparsers.add_parser('deseq_analyse', help='Perform DESeq2 analysis')
    # parser_deseq.add_argument('--gene', required=False, help='Gene name')
    # parser_deseq.set_defaults(func=deseq_analyse)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
