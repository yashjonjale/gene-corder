import argparse
import os
import sys
import json
import subprocess
import requests
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pytximport import tximport
import pyranges as pr
import numpy as np

import pandas as pd
import pyranges as pr
import os
from pyfaidx import Fasta
from Bio import SeqIO



def map_gene_ids_to_names(gtf_file, gene_counts_df):
    """
    Maps gene IDs to gene names using a GTF file and updates a gene counts DataFrame.

    Args:
        gtf_file (str): Path to the GTF file.
        gene_counts_df (pd.DataFrame): A DataFrame with gene IDs in the first column.

    Returns:
        pd.DataFrame: Updated DataFrame with gene names mapped to gene IDs.
    """

    # Step 1: Check if the GTF file exists
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file '{gtf_file}' not found.")
    
    # Step 2: Load the GTF file into a pyranges object
    try:
        gtf = pr.read_gtf(gtf_file)
    except Exception as e:
        raise RuntimeError(f"Error loading GTF file: {e}")
    
    # Step 3: Extract gene_id and gene_name from GTF
    try:
        gtf_genes = gtf[gtf.Feature == 'gene']
        gtf_genes_df = gtf_genes.df[['gene_id', 'gene_name']].drop_duplicates()

        # Step 4: Create a dictionary to map gene_id to gene_name
        gene_id_to_name = dict(zip(gtf_genes_df['gene_id'], gtf_genes_df['gene_name']))

    except KeyError as e:
        raise ValueError(f"Required columns not found in the GTF file: {e}")

    # Step 5: Check if gene_id column exists in gene_counts_df
    if 'gene_id' not in gene_counts_df.columns:
        raise ValueError("The 'gene_id' column is missing from the gene counts dataframe.")
    
    # Step 6: Map gene_id to gene_name in the gene counts DataFrame
    try:
        gene_counts_df['gene_name'] = gene_counts_df['gene_id'].map(gene_id_to_name)
    except Exception as e:
        raise RuntimeError(f"Error mapping gene IDs to gene names: {e}")
    
    # Step 7: Return the updated dataframe
    return gene_counts_df





# 1. Extract Gene Sequence
def extract_gene_sequence(gene_id, gtf_file, genome_fa, output_fa_path):
    # Error handling for file existence
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file '{gtf_file}' not found.")
    if not os.path.exists(genome_fa):
        raise FileNotFoundError(f"Genome FASTA file '{genome_fa}' not found.")
    
    # Load the GTF file
    try:
        gr = pr.read_gtf(gtf_file)
    except Exception as e:
        raise ValueError(f"Error reading GTF file: {str(e)}")
    
    # Filter for gene_name
    # gene_gtf = gr.df[gr.df['gene_name'] == gene_name]
    gene_gtf = gr.df[gr.df['gene_id'] == gene_id]
    if gene_gtf.empty:
        raise ValueError(f"Gene '{gene_id}' not found in the GTF file.")
    
    # Load the genome FASTA
    try:
        genome = Fasta(genome_fa)
    except Exception as e:
        raise ValueError(f"Error loading genome FASTA file: {str(e)}")

    # Initialize sequence
    gene_sequence = ""

    # Extract sequences based on coordinates
    try:
        for _, row in gene_gtf.iterrows():
            chrom = row['Chromosome']
            start = row['Start']
            end = row['End']
            strand = row['Strand']
            
            # Extract sequence
            seq = genome[chrom][start-1:end]  # pyfaidx is 0-based

            if strand == '-':
                seq = seq.reverse.complement

            gene_sequence += str(seq)
    except Exception as e:
        raise ValueError(f"Error extracting sequence: {str(e)}")
    
    # Write to a FASTA file
    try:
        with open(output_fa_path, 'w') as f:
            f.write(f">{gene_id}\n")
            for i in range(0, len(gene_sequence), 60):
                f.write(gene_sequence[i:i+60] + "\n")
        print(f"Gene sequence for {gene_id} saved to {output_fa_path}")
    except Exception as e:
        raise IOError(f"Error writing to FASTA file: {str(e)}")


# 2. Parse GTF for gene mappings
def create_gene_maps(gtf_file):
    # Error handling for file existence
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file '{gtf_file}' not found.")
    
    try:
        gr = pr.read_gtf(gtf_file)
    except Exception as e:
        raise ValueError(f"Error reading GTF file: {str(e)}")

    # Filter and create mappings
    try:
        df = gr.df[['gene_id', 'gene_name']].drop_duplicates()
        gene_id_to_name = dict(zip(df['gene_id'], df['gene_name']))
        gene_name_to_id = dict(zip(df['gene_name'], df['gene_id']))
    except Exception as e:
        raise ValueError(f"Error processing GTF data: {str(e)}")
    
    return gene_id_to_name, gene_name_to_id


# 3. Find gene by sequence
def find_gene_by_sequence(input_fa, gtf_file, genome_fa):
    # Error handling for file existence
    if not os.path.exists(input_fa):
        raise FileNotFoundError(f"Input FASTA file '{input_fa}' not found.")
    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"GTF file '{gtf_file}' not found.")
    if not os.path.exists(genome_fa):
        raise FileNotFoundError(f"Genome FASTA file '{genome_fa}' not found.")
    
    # Load the input sequence
    try:
        input_seq_record = list(SeqIO.parse(input_fa, "fasta"))
        if len(input_seq_record) != 1:
            raise ValueError("Input FASTA file should contain exactly one sequence.")
        input_sequence = str(input_seq_record[0].seq)
    except Exception as e:
        raise ValueError(f"Error reading input FASTA file: {str(e)}")
    
    # Load GTF file
    try:
        gr = pr.read_gtf(gtf_file)
    except Exception as e:
        raise ValueError(f"Error reading GTF file: {str(e)}")

    # Load genome
    try:
        genome = Fasta(genome_fa)
    except Exception as e:
        raise ValueError(f"Error loading genome FASTA file: {str(e)}")
    
    # Dictionary to store gene sequences
    gene_sequences = {}

    # Iterate over each gene in the GTF
    try:
        for gene_name, gene_gtf in gr.df.groupby("gene_name"):
            gene_sequence = ""
            for _, row in gene_gtf.iterrows():
                chrom = row['Chromosome']
                start = row['Start']
                end = row['End']
                strand = row['Strand']

                seq = genome[chrom][start-1:end]  # pyfaidx is 0-based

                if strand == '-':
                    seq = seq.reverse.complement

                gene_sequence += str(seq)

            # Store gene sequence
            gene_sequences[gene_name] = gene_sequence
    except Exception as e:
        raise ValueError(f"Error processing GTF or genome data: {str(e)}")

    # Compare input sequence with each gene's sequence
    try:
        for gene_name, gene_sequence in gene_sequences.items():
            if input_sequence in gene_sequence or gene_sequence in input_sequence:
                print(f"Input sequence matches gene: {gene_name}")
                return gene_name
    except Exception as e:
        raise ValueError(f"Error comparing sequences: {str(e)}")

    print("No matching gene found for the input sequence.")
    return None


def process_gene_counts(obj, abundances_tsv_paths, obj_name, name):
    # Load paths from the object
    t2g_path = obj["index_paths"]["t2g"]
    gtf_path = obj["paths"]["gtf"]
    print(f"[DEBUG] t2g path: {t2g_path}")

    # Load the t2g (transcript-to-gene) map
    t2g_df = pd.read_csv(t2g_path, sep="\t", header=None,
                         names=["transcript_id", "gene_id", "gene_name", "transcript_name", "chrom", "start", "end", "strand"])
    t2g = t2g_df[["transcript_id", "gene_id"]]

    files = abundances_tsv_paths
    print(f"[DEBUG] Files to import: {files}")
    print(f"[DEBUG] Gene map path: {t2g_path}")

    try:
        # Use pytximport to aggregate transcript-level data to gene-level
        result = tximport(files=files, transcript_gene_map=t2g, type='kallisto')
        result = result.round().astype(int)  # Round counts and convert to integer

    except Exception as e:
        print(f"[ERROR] Tximport Failure: {e}")
        return

    # Extract the gene-level counts matrix and gene names
    print(f"[DEBUG] Extracting gene-level counts matrix and gene names")
    
    # Step 1: Load GTF using pyranges
    gtf = pr.read_gtf(gtf_path)
    
    # Step 2: Extract relevant gene_id and gene_name columns from the GTF file
    gtf_gene = gtf.df[gtf.df["Feature"] == "gene"][["gene_id", "gene_name"]].drop_duplicates()
    
    # Step 3: Handle missing gene names by filling them with the gene_id
    gtf_gene["gene_name"].fillna(gtf_gene["gene_id"], inplace=True)
    
    # Step 4: Set gene_id as the index for easy lookups
    gtf_gene.set_index("gene_id", inplace=True)
    
    # Step 5: Find common genes between tximport result and GTF file
    common_genes = result.index.intersection(gtf_gene.index)

    # Step 6: Filter result for common genes
    result = result.loc[common_genes]

    # Step 7: Update gene names based on the GTF file
    result.index = gtf_gene.loc[common_genes, "gene_name"]

    # Step 8: Save the gene-level counts matrix as a CSV file
    print(f"[DEBUG] Saving gene-level counts matrix as a CSV file")
    
    # Create directories if they don't exist
    output_dir = f"./data/{obj_name}/{name}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the gene counts matrix
    gene_count_csv = os.path.join(output_dir, f"{name}_gene_counts.csv")
    result.to_csv(gene_count_csv)

    print(f"[DEBUG] Gene counts saved to {gene_count_csv}")


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

    # constructing the gene maps
    gene_id_to_name, gene_name_to_id = create_gene_maps(gtf_path)

    #save them as json files, and store the paths in the object
    gene_id_to_name_path = os.path.join(output_dir, f"{organism_name}_gene_id_to_name.json")
    gene_name_to_id_path = os.path.join(output_dir, f"{organism_name}_gene_name_to_id.json")

    with open(gene_id_to_name_path, 'w') as f:
        json.dump(gene_id_to_name, f, indent=4)
    with open(gene_name_to_id_path, 'w') as f:
        json.dump(gene_name_to_id, f, indent=4)
    
    obj["gene_maps"] = {
        "gene_id_to_name": gene_id_to_name_path,
        "gene_name_to_id": gene_name_to_id_path
    }


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
    # if name not in obj["quantifications"]:
    print(f"[DEBUG] Creating new quantification entry '{name}'.")
    obj["quantifications"][name] = {
        "sra":{},
        "counts_path": None,
    }
    # else:
    #     print(f"[DEBUG] Quantification entry '{name}' already exists.")
    #     return
    
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
        
        #now download the metadata.tsv for the SRA
        meta_path = os.path.join(sra_dir, f"{sra}/{sra}_metadata.tsv")
        meta_cmd = [
            "pysradb",
            "metadata",
            "--detailed",
            sra
        ]

        print(f"[DEBUG] Running metadata command: {' '.join(meta_cmd)}")
        try:
            # subprocess.run(meta_cmd, check=True)
            with open(meta_path, "w") as meta_file:
                subprocess.run(meta_cmd, check=True, stdout=meta_file)
            print(f"Metadata file {sra} downloaded to {sra_dir}")
            obj["quantifications"][name]["sra"][sra]["meta_path"] = meta_path
        except subprocess.CalledProcessError as e:
            print(f"Error while downloading Metadata file {sra}: {e}")
            return

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
            subprocess.run(fastqdump_cmd, check=True)
            print(f"Fastq files generated for {sra}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running fastq-dump for {sra}: {e}")
            return

    print(f"Fastq dump Done for {sra_list}")

    # ## Now, quantifying the fastq files using kallisto, different commands for paired and single-end reads
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
                "-t", "24",
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
                "-t", "24",
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
        abundances_tsv_path = os.path.join(output_dir, "abundance.tsv")
        abundances_tsv_path = os.path.abspath(abundances_tsv_path)
        abundances_tsv_paths.append(abundances_tsv_path)
        print(f"[DEBUG] Abundance TSV path: {abundances_tsv_path}")


    # update the obj with the new paths
    print("Saving the new paths in config.json\n")
    print(config['objects'][obj_name])
    config['objects'][obj_name] = obj
    save_config(config)


    t2g_path = obj["index_paths"]["t2g"]
    gtf_path = obj["paths"]["gtf"]
    print(f"[DEBUG] t2g path: {t2g_path}")
    
    #load the t2g file as a dataframe
    t2g_df = pd.read_csv(t2g_path, sep="\t", header=None,names=["transcript_id", "gene_id", "gene_name", "transcript_name", "chrom", "start", "end", "strand"])
    t2g = t2g_df[["transcript_id", "gene_id"]]



    files = abundances_tsv_paths
    print(f"[DEBUG] Files to import: {files}")
    print(f"[DEBUG] Gene map path: {t2g_path}")

    try:
        result = tximport(file_paths=files, data_type='kallisto', transcript_gene_map=t2g,
               id_column = "target_id",
               counts_column = "est_counts",
               length_column = "eff_length", abundance_column = "tpm")
        result.X = result.X.round().astype(int)
        # result.obs = pd.DataFrame(index=result.obs_names)
        result.obs["type"] = "type-none"
        result.obs["replicate"] = "replicate-none"

    except Exception as e:
        print(f"[ERROR] Tximport Failure: {e}")
        return

    # Extract the gene-level counts matrix and gene names
    # print(f"[DEBUG] Extracting gene-level counts matrix and gene names")
    # gtf = pr.read_gtf("/data1/yashjonjale/igem_stuff/ehux_study/sra_data/working_folder/kb_ref_out/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf")
    # gtf_gene = gtf.df[gtf.df["Feature"] == "gene"][["gene_id", "gene_name"]].drop_duplicates()
    # gtf_gene["gene_name"].fillna(gtf_gene["gene_id"], inplace=True)
    # gtf_gene.set_index("gene_id", inplace=True)
    # common_genes = result.var.index.intersection(gtf_gene.index)


    # result = result[:, result.var.index.isin(common_genes)]
    # # result.obs = None
    # result.obs = pd.DataFrame(index=result.obs_names)

    # # Save the gene-level counts matrix as a CSV file
    # print(f"[DEBUG] Saving gene-level counts matrix as a CSV file")
    # gene_count_csv = f"./data/{obj_name}/{name}/{name}_gene_counts.csv"
    # result.X = result.X.astype(int)
    # result.var["gene_name"] = gtf_gene.loc[result.var.index, "gene_name"]
    # result.var.set_index("gene_name", inplace=True)
    # result.X = result.X.T
    # result.X.to_csv(gene_count_csv)
    # print(f"[DEBUG] Gene counts saved to {gene_count_csv}")

    # # Update the config with the path to gene count
    # obj["quantifications"][name]["counts_path"] = gene_count_csv

    # print("Saving the new paths in config.json\n")    
    # config['objects'][obj_name] = obj
    # save_config(config)
    # print(f"[DEBUG] Quantification '{name}' completed and saved to config.")



    # print("[DEBUG] quantize function completed successfully.")

    print(f"[DEBUG] Extracting gene-level counts matrix and gene names")
    try:
        # Load GTF and extract gene names
        gtf = pr.read_gtf(gtf_path)
        gtf_gene = gtf.df[gtf.df["Feature"] == "gene"][["gene_id", "gene_name"]].drop_duplicates()
        gtf_gene["gene_name"].fillna(gtf_gene["gene_id"], inplace=True)
        gtf_gene.set_index("gene_id", inplace=True)

        # Find common genes between result and gtf_gene
        common_genes = result.var.index.intersection(gtf_gene.index)

        # Subset to common genes
        result = result[:, result.var.index.isin(common_genes)]

        # Add gene names
        result.var["gene_name"] = gtf_gene.loc[result.var.index, "gene_name"]
        result.var.set_index("gene_name", inplace=True)

        file_to_sra_mapping = dict(zip(abundances_tsv_paths, sra_list))

        # Replace result.obs_names with the SRA codes
        result.obs_names = result.obs_names.map(lambda x: file_to_sra_mapping.get(x, x))  # Replace path with SRA code

        # Print the updated result.obs_names for debugging
        print(f"[DEBUG] result.obs_names (after replacing with SRA codes): {result.obs_names}")

    except Exception as e:
        print(f"[ERROR] Error in extracting gene names: {e}")
        return

    # # Step 4: Save gene-level counts matrix to CSV
    # try:
    #     print(f"[DEBUG] Saving gene-level counts matrix as a CSV file")
    #     gene_count_csv = os.path.join(f"./data/{obj_name}/{name}", f"{name}_gene_counts.csv")

    #     os.makedirs(os.path.dirname(gene_count_csv), exist_ok=True)

    #     result.X = result.X.astype(int).T  # Transpose to make genes as columns
    #     result.X.to_csv(gene_count_csv)

    #     print(f"[DEBUG] Gene counts saved to {gene_count_csv}")

    #     # Update config with the path to the gene counts
    #     obj["quantifications"][name]["counts_path"] = gene_count_csv

    #     print(f"Saving the new paths in config.json")
    #     config['objects'][obj_name] = obj
    #     save_config(config)
    #     print(f"[DEBUG] Quantification '{name}' completed and saved to config.")

    # except Exception as e:
    #     print(f"[ERROR] Error while saving gene-level counts matrix: {e}")
    #     return


    try:
        print(f"[DEBUG] Saving gene-level counts matrix as a CSV file")
        gene_count_csv = os.path.join(f"./data/{obj_name}/{name}", f"{name}_gene_counts.csv")
        os.makedirs(os.path.dirname(gene_count_csv), exist_ok=True)

        # Convert the numpy array result.X to a pandas DataFrame
        gene_counts_df = pd.DataFrame(result.X, index=result.obs_names, columns=result.var_names)

        print(f"[DEBUG] result,obs_names: {result.obs_names}")

        print(f"[DEBUG] Gene counts DataFrame shape: {gene_counts_df.shape}")
        print(f"[DEBUG] Gene counts DataFrame head:\n{gene_counts_df.head()}")
        #columns names are 
        print(f"[DEBUG] Gene counts DataFrame columns:\n{gene_counts_df.columns}")
        # Transpose to make genes as columns
        gngdf = gene_counts_df.T
        print(f"[DEBUG] Gene counts DataFrame shape after transposing: {gngdf.shape}")
        print(f"[DEBUG] Gene counts DataFrame head after transposing:\n{gngdf.head()}")
        print(f"[DEBUG] Gene counts DataFrame columns after transposing:\n{gngdf.columns}") 
        # Save the gene counts matrix to a CSV file
        gngdf.to_csv(gene_count_csv)

        print(f"[DEBUG] Gene counts saved to {gene_count_csv}")

        # Update config with the path to the gene counts
        obj["quantifications"][name]["counts_path"] = gene_count_csv

        print(f"Saving the new paths in config.json")
        config['objects'][obj_name] = obj
        save_config(config)
        print(f"[DEBUG] Quantification '{name}' completed and saved to config.")

    except Exception as e:
        print(f"[ERROR] Error while saving gene-level counts matrix: {e}")
        return
    print("[DEBUG] quantize function completed successfully.")

def list_quantifications(args):
    obj_name = args.obj
    print(f"Listing quantifications for object '{obj_name}'...")
    config = None
    print("[DEBUG] Loading config.json for listing quantifications.")
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

    if not obj["quantifications"]:
        print(f"[DEBUG] No quantifications found for object '{obj_name}'.")
        return
    
    print(f"Quantifications for object '{obj_name}':")
    for quant_name, quant in obj["quantifications"].items():
        print(f"- {quant_name}:")
        for sra, sra_data in quant["sra"].items():
            print(f"  - {sra}:")
            print(f"    - Path: {sra_data['path']}")
            print(f"    - Paired: {sra_data['paired']}")
            print(f"    - Fastq Dir: {sra_data['fastq_dir']}")
            print(f"    - Quant Path: {sra_data['quant_path']}")
            print(f"    - Metadata Path: {sra_data['meta_path']}")
        print(f"  - Counts Path: {quant['counts_path']}")
    print("[DEBUG] list_quantifications function completed successfully.")


def plot_gene_abundances(args):
    obj_name = args.obj
    quant_name = args.quantification_name
    target_gene = args.gene
    named = args.named

    if named:
        #find the gene_id from the gene_name
        gene_id = None
        with open(f"./data/{obj_name}/{quant_name}/{obj_name}_gene_name_to_id.json") as f:
            gene_name_to_id = json.load(f)
            gene_id = gene_name_to_id.get(target_gene, None)
        if gene_id is None:
            print(f"Gene name '{target_gene}' not found in the gene name to ID mapping.")
            return
        else:
            target_gene = gene_id
    # Debug output
    print(f"[DEBUG] Object Name: {obj_name}")
    config_path = 'config.json'

    # Step 1: Load the configuration
    try:
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file '{config_path}' not found.")

        with open(config_path, 'r') as f:
            config = json.load(f)

        if obj_name not in config['objects']:
            raise ValueError(f"Object '{obj_name}' not found in config.")
            
        obj = config['objects'][obj_name]
        
        if quant_name not in obj["quantifications"]:
            raise ValueError(f"Quantification '{quant_name}' not found for object '{obj_name}'.")
            
        # Step 2: Get the path to the gene counts file
        path_to_gene_count = obj["quantifications"][quant_name]["counts_path"]
        if not os.path.exists(path_to_gene_count):
            raise FileNotFoundError(f"Gene counts file '{path_to_gene_count}' not found.")
        
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"[ERROR] {e}")
        return

    # Step 3: Load gene counts
    try:
        print(f"[DEBUG] Loading gene counts from '{path_to_gene_count}'")
        df = pd.read_csv(path_to_gene_count)
        if df.empty:
            raise ValueError(f"The gene counts file '{path_to_gene_count}' is empty.")
        print(f"OK")    
        # Ensure the target gene exists in the dataframe
        if target_gene not in df.iloc[:, 0].values:
            raise ValueError(f"Target gene '{target_gene}' not found in gene counts file.")
        
    except pd.errors.EmptyDataError:
        print(f"[ERROR] The gene counts file '{path_to_gene_count}' is empty or corrupt.")
        return
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        return
    except Exception as e:
        print(f"[ERROR] An error occurred while loading the gene counts: {e}")
        return

    # Step 4: Extract the gene row and plot data
    try:
        # Convert the dataframe to numpy for fast lookup
        array_with_genes = df.values
        
        # Fetch the row where the first column matches the target gene
        gene_row = array_with_genes[array_with_genes[:, 0] == target_gene]
        
        if gene_row.size == 0:
            raise ValueError(f"Gene '{target_gene}' not found in the gene counts file.")
        
        # Extract the gene's expression values across samples
        values = gene_row[0, 1:].astype(float)  # Assuming gene expression values are numeric
        labels = df.columns[1:]  # Sample names are in the columns
        
        # Step 5: Plot the gene abundances
        plt.figure(figsize=(10, 6))
        plt.bar(labels, values, color='blue')
        plt.xlabel("Samples")
        plt.ylabel("Expression Level (in TPM)")
        plt.title(f"Gene Abundances for {target_gene}")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig("output.png")

    except IndexError as e:
        print(f"[ERROR] Index error occurred while processing gene row: {e}")
    except ValueError as e:
        print(f"[ERROR] Value error: {e}")
    except Exception as e:
        print(f"[ERROR] An unexpected error occurred during plotting: {e}")


def generate_correlation_matrix(args):
    try:
        # Extract arguments
        obj_name = args.obj
        quant_name = args.quantification_name
        gene_list_file = args.genes  # Path to gene list file, txt file with gene ids

        # Debug output
        print(f"Object Name: {obj_name}, Quantification Name: {quant_name}, Gene List File: {gene_list_file}")

        # Load the config file
        if not os.path.exists('config.json'):
            raise FileNotFoundError("The config.json file was not found.")

        with open('config.json', 'r') as f:
            try:
                config = json.load(f)
            except json.JSONDecodeError as e:
                raise ValueError(f"Error decoding config.json: {str(e)}")

        # Retrieve object and path to gene count
        if obj_name not in config.get('objects', {}):
            raise KeyError(f"Object '{obj_name}' not found in config.json.")

        obj = config['objects'][obj_name]
        quantifications = obj.get('quantifications', {})

        if quant_name not in quantifications:
            raise KeyError(f"Quantification '{quant_name}' not found for object '{obj_name}' in config.json.")

        path_to_gene_count = quantifications[quant_name].get("counts_path")
        if not os.path.exists(path_to_gene_count):
            raise FileNotFoundError(f"Gene count file '{path_to_gene_count}' not found.")

        # Load gene count data
        try:
            gene_count_df = pd.read_csv(path_to_gene_count)
        except Exception as e:
            raise ValueError(f"Error reading gene count file '{path_to_gene_count}': {str(e)}")



        gene_count_data = gene_count_df.values
        
        #debug output
        print(f"[DEBUG] Gene count data shape: {gene_count_data.shape}")
        print(f"[DEBUG] Gene count data head:\n{gene_count_data[:5]}")

        # Load gene list file
        if not os.path.exists(gene_list_file):
            raise FileNotFoundError(f"Gene list file '{gene_list_file}' not found.")
        
        try:
            gene_list_df = pd.read_csv(gene_list_file)
        except Exception as e:
            raise ValueError(f"Error reading gene list file '{gene_list_file}': {str(e)}")

        gene_list = gene_list_df.values

        # Filter gene count data based on the gene list
        try:
            array_with_genes = gene_count_data[np.isin(gene_count_data[:, 0], gene_list[:, 0])]
            if array_with_genes.size == 0:
                raise ValueError("No matching genes found between the gene count data and the provided gene list.")
        except Exception as e:
            raise ValueError(f"Error filtering gene count data: {str(e)}")

        # Extract gene names (first column)
        gene_names = array_with_genes[:, 0]

        # Extract the numerical data (excluding the first column)
        try:
            data = array_with_genes[:, 1:].astype(float)
        except ValueError:
            raise ValueError("Non-numeric data found in the gene count file, unable to convert to float.")

        # Create a DataFrame with gene names as index
        df = pd.DataFrame(data, index=gene_names)

        # Compute the correlation matrix
        correlation_matrix = df.corr()

        # Plot the heatmap
        plt.figure(figsize=(10, 8))  # Adjust figure size for better readability
        sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', 
                    xticklabels=gene_names, yticklabels=gene_names, 
                    annot_kws={"size": 8},  # Adjust font size
                    cbar_kws={'label': 'Correlation coefficient'})  # Add colorbar label

        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right', fontsize=10)
        plt.yticks(fontsize=10)

        # Adjust layout
        plt.tight_layout()

        # Display the plot
        plt.show()

    except FileNotFoundError as fnf_error:
        print(f"File error: {str(fnf_error)}")
    except KeyError as key_error:
        print(f"Key error: {str(key_error)}")
    except ValueError as value_error:
        print(f"Value error: {str(value_error)}")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")


def name2id(args):
    obj_name = args.obj
    gene_name = args.gene_name
    print(f"Converting gene name '{gene_name}' to gene ID for object '{obj_name}'...")
    config = None
    try:
        with open('config.json', 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print("[ERROR] config.json not found. Please instantiate an organism first.")
        return

    if obj_name not in config['objects']:
        print(f"[ERROR] Object '{obj_name}' not found in config.")
        return

    obj = config['objects'][obj_name]
    print(f"[DEBUG] Loaded object '{obj_name}' from config.")

    gene_id = None
    gene_name_to_id_path = obj["gene_maps"]["gene_name_to_id"]
    if not os.path.exists(gene_name_to_id_path):
        print(f"[ERROR] Gene name to ID mapping file '{gene_name_to_id_path}' not found.")
        return

    try:
        with open(gene_name_to_id_path, 'r') as f:
            gene_name_to_id = json.load(f)
            gene_id = gene_name_to_id.get(gene_name, None)
    except Exception as e:
        print(f"[ERROR] Error loading gene name to ID mapping: {e}")
        return

    if gene_id is None:
        print(f"Gene name '{gene_name}' not found in the gene name to ID mapping.")
    else:
        print(f"Gene name '{gene_name}' maps to gene ID '{gene_id}'.")

def gene2fasta(args):
    obj_name = args.obj
    gene_name = args.gene
    named = args.named
    gene_id = None
    output_dir = args.output_dir
    # print(f"Converting gene name '{gene_name}' to gene ID for object '{obj_name}'...")
    config = None
    try:
        with open('config.json', 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print("[ERROR] config.json not found. Please instantiate an organism first.")
        return

    if obj_name not in config['objects']:
        print(f"[ERROR] Object '{obj_name}' not found in config.")
        return

    obj = config['objects'][obj_name]
    
    print(f"[DEBUG] Loaded object '{obj_name}' from config.")

    if named:
        gene_name_to_id_path = obj["gene_maps"]["gene_name_to_id"]
        if not os.path.exists(gene_name_to_id_path):
            print(f"[ERROR] Gene name to ID mapping file '{gene_name_to_id_path}' not found.")
            return

        try:
            with open(gene_name_to_id_path, 'r') as f:
                gene_name_to_id = json.load(f)
                gene_id = gene_name_to_id.get(gene_name, None)
        except Exception as e:
            print(f"[ERROR] Error loading gene name to ID mapping: {e}")
            return

        if gene_id is None:
            print(f"Gene name '{gene_name}' not found in the gene name to ID mapping.")
        else:
            print(f"Gene name '{gene_name}' maps to gene ID '{gene_id}'.")
        
    else:
        gene_id = gene_name
    
    outpath = os.path.join(output_dir, f"{gene_name}.fasta")
    gtf_path = obj["paths"]["gtf"]
    genome_path = obj["paths"]["genome"]    

    # Ensure GTF and genome paths exist
    if not os.path.exists(gtf_path):
        print(f"[ERROR] GTF file '{gtf_path}' not found.")
        return
    if not os.path.exists(genome_path):
        print(f"[ERROR] Genome file '{genome_path}' not found.")
        return

    # At this point, all inputs are valid, and you can call the function to extract the gene sequence
    try:
        extract_gene_sequence(gene_id, gtf_path, genome_path, outpath)
    except Exception as e:
        print(f"[ERROR] Failed to extract gene sequence: {e}")
    


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
    parser_list_quant = subparsers.add_parser('list_quant', help='List quantifications')
    parser_list_quant.add_argument('--obj', required=True, help='Object name')
    parser_list_quant.set_defaults(func=list_quantifications)

    parser_plot = subparsers.add_parser('plot_gene_abundances', help='Plot gene abundances')
    parser_plot.add_argument('--gene', required=False, help='Gene name')
    parser_plot.add_argument('--named', required=False, action='store_true', help='Gene name is provided')
    parser_plot.add_argument('--obj', required=True, help='Object name')
    parser_plot.add_argument('--quantification_name', required=True, help='Quantification name')
    parser_plot.set_defaults(func=plot_gene_abundances)

    parser_corr = subparsers.add_parser('generate_correlation_matrix', help='Generate correlation matrix')
    parser_corr.add_argument('--genes', required=False, help='newline separated list of genes')
    parser_corr.add_argument('--obj', required=True, help='Object name')
    parser_corr.add_argument('--quantification_name', required=True, help='Quantification name')
    parser_corr.set_defaults(func=generate_correlation_matrix)

    parser_name2id = subparsers.add_parser('name2id', help='Convert gene name to gene ID')
    parser_name2id.add_argument('--gene_name', required=True, help='Gene name')
    parser_name2id.add_argument('--obj', required=True, help='Object name')
    parser_name2id.set_defaults(func=name2id)

    parser_gene2fasta = subparsers.add_parser('gene2fasta', help='Extract gene sequence to FASTA')
    parser_gene2fasta.add_argument('--gene', required=True, help='Gene name or ID')
    parser_gene2fasta.add_argument('--obj', required=True, help='Object name')
    parser_gene2fasta.add_argument('--output_dir', required=True, help='Output directory')
    parser_gene2fasta.add_argument('--named', required=False, action='store_true', help='Gene name is provided')

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
