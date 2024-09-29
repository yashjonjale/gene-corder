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
    ## Save the config file
    print("[DEBUG] Saving configuration to config.json")
    name = 'config.json'
    # Delete the existing file config.json
    if os.path.exists(name):
        print(f"[DEBUG] Existing {name} found. Removing it.")
        os.remove(name)
    else:
        print(f"[DEBUG] No existing {name} found. Creating a new one.")
    with open(name, 'w') as f:
        json.dump(config_json, f, indent=4)
    print("[DEBUG] Configuration saved successfully.")

def instantiate_organism(args):
    print("[DEBUG] Starting instantiate_organism function.")
    print(f"[DEBUG] Received arguments: {args}")

    print(f"Instantiating Organism ....")

    organism_name = args.organism
    name = args.name
    desc = args.desc
    transcript_path = args.transcriptome_path
    genome_path = args.genome_path
    gtf_path = args.gtf_path

    print(f"Organism: {organism_name}")
    print(f"Name: {name}")
    print(f"Instantiating object '{name}' for organism '{organism_name}'...")
    # Debug output
    print(f"Transcriptome Path: {transcript_path}")
    print(f"Genome Path: {genome_path}")
    print(f"GTF Path: {gtf_path}")

    # Check if paths are absolute, if not, make them absolute
    if transcript_path is not None:
        if not os.path.isabs(transcript_path):
            transcript_path = os.path.abspath(transcript_path)
            print(f"[DEBUG] Converted transcriptome_path to absolute path: {transcript_path}")
    if not os.path.isabs(genome_path):
        genome_path = os.path.abspath(genome_path)
        print(f"[DEBUG] Converted genome_path to absolute path: {genome_path}")

    if not os.path.isabs(gtf_path):
        gtf_path = os.path.abspath(gtf_path)
        print(f"[DEBUG] Converted gtf_path to absolute path: {gtf_path}")

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
            # Placeholder for quantifications
        },
        "index_paths": {
            "kallisto_index": None,
            "t2g": None
        }
    }

    # Load the current config and append the new object
    print("[DEBUG] Loading existing config.json")
    config = None
    try:
        with open('config.json', 'r') as f:
            config = json.load(f)
            print(f"[DEBUG] Current config objects: {list(config['objects'].keys())}")
    except FileNotFoundError:
        print("[DEBUG] config.json not found. Initializing a new config.")
        config = {"objects": {}}

    config['objects'][name] = obj
    print(f"[DEBUG] Added object '{name}' to config.")

    # Directory to store outputs
    output_dir = f"./data/{name}"
    print(f"[DEBUG] Creating output directory at {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    ## Run gffread to convert GTF to FASTA if transcript_path is not provided or doesn't exist
    if transcript_path is None or not os.path.exists(transcript_path):
        print(f"Converting GTF to FASTA...")
        gffread_cmd = [
            "gffread", "-w", f"{output_dir}/{name}_transcripts.fa", "-g", genome_path, gtf_path
        ]
        print(f"[DEBUG] Running gffread command: {' '.join(gffread_cmd)}")
        try:
            subprocess.run(gffread_cmd, check=True)
            print(f"FASTA file generated at {output_dir}/{name}_transcripts.fa")
            obj["paths"]["transcriptome"] = f"{output_dir}/{name}_transcripts.fa"
            transcript_path = f"{output_dir}/{name}_transcripts.fa"
            # Make absolute path
            transcript_path = os.path.abspath(transcript_path)
            print(f"Transcript FastA path: {transcript_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running gffread: {e}")
            return
    else:
        print(f"[DEBUG] Transcriptome path exists: {transcript_path}")

    # Step 1: Index genome using kallisto index
    kallisto_index_path = os.path.join(output_dir, f"index.idx")
    kallisto_index_path = os.path.abspath(kallisto_index_path)
    print(f"Indexing genome with kallisto...")
    kallisto_cmd = [
        "kallisto", "index",
        "-i", kallisto_index_path,
        transcript_path
    ]
    print(f"[DEBUG] Running kallisto index command: {' '.join(kallisto_cmd)}")
    
    try:
        subprocess.run(kallisto_cmd, check=True)
        print(f"Kallisto index created at {kallisto_index_path}")
        obj["index_paths"]["kallisto_index"] = kallisto_index_path
    except subprocess.CalledProcessError as e:
        print(f"Error while running kallisto index: {e}")
        return

    # Step 2: Generate t2g.tsv using kb ref
    t2g_path = os.path.join(output_dir, f"{name}_t2g.tsv")
    t2g_path = os.path.abspath(t2g_path)
    print(f"Generating t2g.tsv using kb ref...")
    kb_ref_cmd = [
        "kb", "ref", "-i", "ind.idx",
        "-g", t2g_path,
        "-f1", transcript_path,
        genome_path, gtf_path
    ]
    print(f"[DEBUG] Running kb ref command: {' '.join(kb_ref_cmd)}")
    
    try:
        subprocess.run(kb_ref_cmd, check=True)
        print(f"t2g.tsv generated at {t2g_path}")
        obj["index_paths"]["t2g"] = t2g_path
        # Delete ind.idx
        if os.path.exists("ind.idx"):
            print("[DEBUG] Deleting ind.idx file.")
            os.remove("ind.idx")
        else:
            print("[DEBUG] ind.idx file does not exist. Skipping deletion.")
    except subprocess.CalledProcessError as e:
        print(f"Error while running kb ref: {e}")
        return

    # Step 3: Save the updated config file with the new paths
    print("[DEBUG] Updating config object with new paths.")
    print("Updated object:")
    print(json.dumps(obj, indent=4))
    config['objects'][name] = obj
    save_config(config)
    print(f"Updated config file with new paths for '{name}'.")
    
    print(f"Genome indexing and t2g generation complete. Ready for downstream quantification!")
    print("[DEBUG] instantiate_organism function completed successfully.")

def quantize(args):
    print("[DEBUG] Starting quantize function.")
    print(f"[DEBUG] Received arguments: {args}")

    sra_list = args.sra.split(',') if args.sra else []
    obj_name = args.obj
    name = args.name
    paired = args.paired

    print("Quantizing RNA-Seq data...")
    
    # Debug output
    print(f"Name: {name}")
    print(f"SRA List: {sra_list}")
    print(f"Object Name: {obj_name}")
    print(f"Paired-end reads: {paired}")

    config = None
    print("[DEBUG] Loading config.json for quantization.")
    try:
        with open('config.json', 'r') as f:
            config = json.load(f)
            print(f"[DEBUG] Available objects in config: {list(config['objects'].keys())}")
    except FileNotFoundError:
        print("[ERROR] config.json not found. Please instantiate an organism first.")
        return

    if obj_name not in config['objects']:
        print(f"[ERROR] Object '{obj_name}' not found in config.")
        return

    obj = config['objects'][obj_name]
    print(f"[DEBUG] Loaded object '{obj_name}' from config.")

    ## Initialize quantification entry if not present
    if name not in obj["quantifications"]:
        print(f"[DEBUG] Creating new quantification entry '{name}'.")
        obj["quantifications"][name] = {
            "sra":{},
            "counts_path": None,
        }
    else:
        print(f"[DEBUG] Quantification entry '{name}' already exists.")
        return
    
    sra_dir = f"./data/{obj_name}/{name}/"
    os.makedirs(sra_dir, exist_ok=True)
    ## Download SRA files
    for sra in sra_list:
        print(f"[DEBUG] Processing SRA accession: {sra}")
        # Create path for the SRA file
        obj["quantifications"][name]["sra"][sra] = {
            "path": os.path.join(sra_dir, f"{sra}/{sra}.sra"),
            "fastq_dir": os.path.join(sra_dir, f"{sra}/"),
            "paired": paired
        }
        sra_path = os.path.join(sra_dir, f"{sra}/{sra}.sra")
        print(f"[DEBUG] SRA file path: {sra_path}")

        download = False  # Flag to determine if download is required

        # Check if SRA file already exists
        if not os.path.exists(sra_path):
            download = True
            print(f"[DEBUG] SRA file {sra_path} does not exist and will be downloaded.")
        else:
            print(f"[DEBUG] SRA file {sra_path} already exists. Skipping download.")

        # Download the SRA file using prefetch if required
        if download:
            print(f"Downloading {sra}...")
            prefetch_cmd = [
                "prefetch",
                sra,
                "-O",
                sra_dir
            ]
            print(f"[DEBUG] Running prefetch command: {' '.join(prefetch_cmd)}")
            try:
                subprocess.run(prefetch_cmd, check=True)
                print(f"SRA file {sra} downloaded to {sra_dir}")
            except subprocess.CalledProcessError as e:
                print(f"Error while downloading SRA file {sra}: {e}")
                return
        else:
            print(f"[DEBUG] Skipping download for {sra} as it already exists.")

    print(f"SRA prefetch Done for {sra_list}")

    ## Once downloaded, perform fastq-dump on the corresponding SRA files
    for sra in sra_list:
        print(f"[DEBUG] Running fastq-dump for SRA accession: {sra}")
        sra_path = f"./data/{obj_name}/{name}/{sra}/{sra}.sra"
        fastqdump_path = f"./data/{obj_name}/{name}/{sra}/"
        print(f"[DEBUG] SRA path: {sra_path}")
        print(f"[DEBUG] Fastq output directory: {fastqdump_path}")

        if paired:
            fastqdump_cmd = [
                "fastq-dump",
                "--split-files",
                "-O",
                fastqdump_path,
                sra_path
            ]
            print(f"[DEBUG] Running paired-end fastq-dump command: {' '.join(fastqdump_cmd)}")
        else:
            fastqdump_cmd = [
                "fastq-dump",
                "-O",
                fastqdump_path,
                sra_path
            ]
            print(f"[DEBUG] Running single-end fastq-dump command: {' '.join(fastqdump_cmd)}")
        
        try:
            # Uncomment the next line to enable fastq-dump execution
            # subprocess.run(fastqdump_cmd, check=True)
            print(f"Fastq files generated for {sra}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running fastq-dump for {sra}: {e}")
            return

    print(f"Fastq dump Done for {sra_list}")

    ## Now, quantifying the fastq files using kallisto, different commands for paired and single-end reads
    print("Quantifying the fastq files using kallisto...")
    
    # Run kallisto quantification for each fastq file
    abundances_tsv_paths = []
    for sra in sra_list:
        print(f"[DEBUG] Quantifying SRA accession: {sra}")
        fastqdump_path = f"./data/{obj_name}/{name}/{sra}/"
        kallisto_index_path = obj["index_paths"]["kallisto_index"]
        output_dir = f"./data/{obj_name}/{name}/{sra}_quant/"
        os.makedirs(output_dir, exist_ok=True)
        print(f"[DEBUG] Kallisto index path: {kallisto_index_path}")
        print(f"[DEBUG] Kallisto output directory: {output_dir}")
        obj["quantifications"][name]["sra"][sra]["quant_path"] = output_dir
        if paired:
            kallisto_cmd = [
                "kallisto", "quant",
                "-i", kallisto_index_path,
                "-o", output_dir,
                f"{fastqdump_path}/{sra}_1.fastq",
                f"{fastqdump_path}/{sra}_2.fastq"
            ]
            print(f"[DEBUG] Running paired-end kallisto quant command: {' '.join(kallisto_cmd)}")
        else:
            kallisto_cmd = [
                "kallisto", "quant",
                "-i", kallisto_index_path,
                "--single",
                "-l", "200",
                "-s", "20",
                "-o", output_dir,
                f"{fastqdump_path}/{sra}.fastq"
            ]
            print(f"[DEBUG] Running single-end kallisto quant command: {' '.join(kallisto_cmd)}")
        
        try:
            print(f"Running kallisto quant for {sra}...")
            print(kallisto_cmd)
            subprocess.run(kallisto_cmd, check=True)
            print(f"Quantification completed for {sra}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running kallisto quant for {sra}: {e}")
            return

        ## Save the path to the quantification in the obj quantification path entries
        print(f"[DEBUG] Saving quantification results for {sra}")
        obj["quantifications"][name]["paths"].append(os.path.abspath(output_dir))
        obj["quantifications"][name]["types"].append("paired" if paired else "single")
        abundances_tsv_path = os.path.join(output_dir, "abundance.tsv")
        abundances_tsv_path = os.path.abspath(abundances_tsv_path)
        abundances_tsv_paths.append(abundances_tsv_path)
        print(f"[DEBUG] Abundance TSV path: {abundances_tsv_path}")


    # update the obj with the new paths
    print("Saving the new paths in config.json\n")
    print(config['objects'][obj_name])
    config['objects'][obj_name] = obj
    save_config(config)



    # Importing and aggregating gene counts using tximport
    print("[DEBUG] Importing and aggregating gene counts using tximport.")
    try:
        import txi  # Ensure that the 'txi' module is available
    except ImportError:
        print("[ERROR] The 'txi' module is not installed or not found.")
        return

    files = abundances_tsv_paths
    gene_map = obj["index_paths"]["t2g"]  # A file mapping transcripts to genes

    print(f"[DEBUG] Files to import: {files}")
    print(f"[DEBUG] Gene map path: {gene_map}")

    try:
        # Load the transcript data
        txi_data = txi.load_transcript_data(files, gene_map)
        print("[DEBUG] Transcript data loaded successfully.")
    except Exception as e:
        print(f"[ERROR] Failed to load transcript data: {e}")
        return

    # Extract the gene-level counts matrix and gene names
    gene_counts = txi_data.get('counts')
    gene_names = txi_data.get('gene_names')

    if gene_counts is None or gene_names is None:
        print("[ERROR] 'counts' or 'gene_names' not found in txi_data.")
        return

    print("[DEBUG] Gene counts and gene names extracted.")

    # Convert the counts matrix to a DataFrame
    print("[DEBUG] Converting gene counts to DataFrame.")
    df_gene_counts = pd.DataFrame(gene_counts, index=gene_names, columns=sra_list)
    print("[DEBUG] DataFrame created successfully.")

    # Reset the index to make gene names a column instead of the index
    df_gene_counts.reset_index(inplace=True)

    # Rename the first column as 'Gene'
    df_gene_counts.rename(columns={'index': 'Gene'}, inplace=True)
    print("[DEBUG] Renamed the first column to 'Gene'.")

    # Define the gene count directory
    gene_count_dir = f"./data/{obj_name}/{name}/gene_count/"
    os.makedirs(gene_count_dir, exist_ok=True)
    print(f"[DEBUG] Created gene count directory at {gene_count_dir}")

    # Save the DataFrame as a CSV file
    gene_count_csv = os.path.join(gene_count_dir, "gene_count.csv")
    df_gene_counts.to_csv(gene_count_csv, index=False)
    print(f"[DEBUG] Gene counts saved to {gene_count_csv}")

    # Update the config with the path to gene count
    obj["quantifications"][name]["path_to_gene_count"] = gene_count_csv

    print("Saving the new paths in config.json\n")    
    config['objects'][obj_name] = obj
    save_config(config)
    print(f"[DEBUG] Quantification '{name}' completed and saved to config.")

    print("[DEBUG] quantize function completed successfully.")

def main():

    # Check if config file exists
    if not os.path.exists('config.json'):
        print("[DEBUG] config.json does not exist. Creating a new one.")
        with open('config.json', 'w') as f:
            json.dump({"objects": {}}, f, indent=4)
    else:
        print("[DEBUG] config.json found.")

    parser = argparse.ArgumentParser(description='GENECorder CLI')
    subparsers = parser.add_subparsers(dest='command')

    # Instantiate
    parser_instantiate = subparsers.add_parser('instantiate', help='Instantiate an object for an organism')
    parser_instantiate.add_argument('--organism', required=True, help='Organism name')
    parser_instantiate.add_argument('--name', required=True, help='Name for this quantification')
    parser_instantiate.add_argument('--desc', required=False, help='Description for this quantification')
    parser_instantiate.add_argument('--transcriptome_path', required=False, help='Path to transcriptome file in fasta format')
    parser_instantiate.add_argument('--genome_path', required=True, help='Path to genome file in fasta format')
    parser_instantiate.add_argument('--gtf_path', required=True, help='Path to the annotation file in GTF/GFF format')
    parser_instantiate.set_defaults(func=instantiate_organism)

    # Quantize
    parser_quantize = subparsers.add_parser('quantize', help='Quantize RNA-Seq data')
    parser_quantize.add_argument('--sra', required=True, help='Comma-separated SRA accession codes')
    parser_quantize.add_argument('--name', required=True, help='Name for this quantification')
    parser_quantize.add_argument('--obj', required=True, help='Object name')
    parser_quantize.add_argument('--paired', required=False, action='store_true', help='Paired-end reads')
    parser_quantize.set_defaults(func=quantize)

    # Uncomment and implement additional subparsers as needed
    # parser_list_quant = subparsers.add_parser('list_quant', help='List quantifications')
    # parser_list_quant.set_defaults(func=list_quantifications)

    # parser_plot = subparsers.add_parser('plot_gene_abundances', help='Plot gene abundances')
    # parser_plot.add_argument('--gene', required=False, help='Gene name')
    # parser_plot.set_defaults(func=plot_gene_abundances)

    # parser_corr = subparsers.add_parser('generate_correlation_matrix', help='Generate correlation matrix')
    # parser_corr.add_argument('--genes', required=False, help='Comma-separated list of genes')
    # parser_corr.set_defaults(func=generate_correlation_matrix)

    # parser_deseq = subparsers.add_parser('deseq_analyse', help='Perform DESeq2 analysis')
    # parser_deseq.add_argument('--gene', required=False, help='Gene name')
    # parser_deseq.set_defaults(func=deseq_analyse)

    args = parser.parse_args()
    print(f"[DEBUG] Parsed arguments: {args}")

    if hasattr(args, 'func'):
        print(f"[DEBUG] Executing command: {args.command}")
        args.func(args)
    else:
        print("[DEBUG] No command provided. Displaying help.")
        parser.print_help()

if __name__ == "__main__":
    main()
