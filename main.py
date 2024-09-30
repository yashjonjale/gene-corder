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


import pandas as pd
import pyranges as pr
import os

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

# # Example usage:

# gtf_file = "your_file.gtf"  # Replace with the actual path to your GTF file

# # Sample RNA-seq gene counts dataframe
# data = {
#     "gene_id": ["EMIHUDRAFT_100000", "EMIHUDRAFT_100004", "EMIHUDRAFT_100007", "EMIHUDRAFT_100011", "EMIHUDRAFT_100015"],
#     "sample1": [246, 457, 152, 40.5625, 60.7854]
# }
# df = pd.DataFrame(data)

# Call the function and handle potential errors
# try:
#     updated_df = map_gene_ids_to_names(gtf_file, df)
#     print(updated_df)
# except (FileNotFoundError, ValueError, RuntimeError) as e:
#     print(f"[ERROR] {e}")





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
    # print(f"[DEBUG] Creating new quantification entry '{name}'.")
    # obj["quantifications"][name] = {
    #     "sra":{},
    #     "counts_path": None,
    # }
    # # else:
    # #     print(f"[DEBUG] Quantification entry '{name}' already exists.")
    # #     return
    
    # sra_dir = f"./data/{obj_name}/{name}/"
    # os.makedirs(sra_dir, exist_ok=True)
    # ## Download SRA files
    # for sra in sra_list:
    #     print(f"[DEBUG] Processing SRA accession: {sra}")
    #     # Create path for the SRA file
    #     obj["quantifications"][name]["sra"][sra] = {
    #         "path": os.path.join(sra_dir, f"{sra}/{sra}.sra"),
    #         "fastq_dir": os.path.join(sra_dir, f"{sra}/"),
    #         "paired": paired
    #     }
    #     sra_path = os.path.join(sra_dir, f"{sra}/{sra}.sra")
    #     print(f"[DEBUG] SRA file path: {sra_path}")

    #     download = False  # Flag to determine if download is required

    #     # Check if SRA file already exists
    #     if not os.path.exists(sra_path):
    #         download = True
    #         print(f"[DEBUG] SRA file {sra_path} does not exist and will be downloaded.")
    #     else:
    #         print(f"[DEBUG] SRA file {sra_path} already exists. Skipping download.")

    #     # Download the SRA file using prefetch if required
    #     if download:
    #         print(f"Downloading {sra}...")
    #         prefetch_cmd = [
    #             "prefetch",
    #             sra,
    #             "-O",
    #             sra_dir
    #         ]
    #         print(f"[DEBUG] Running prefetch command: {' '.join(prefetch_cmd)}")
    #         try:
    #             subprocess.run(prefetch_cmd, check=True)
    #             print(f"SRA file {sra} downloaded to {sra_dir}")
    #         except subprocess.CalledProcessError as e:
    #             print(f"Error while downloading SRA file {sra}: {e}")
    #             return
    #     else:
    #         print(f"[DEBUG] Skipping download for {sra} as it already exists.")
        
    #     #now download the metadata.tsv for the SRA
    #     meta_path = os.path.join(sra_dir, f"{sra}/{sra}_metadata.tsv")
    #     meta_cmd = [
    #         "pysradb",
    #         "metadata",
    #         "--detailed",
    #         sra
    #     ]

    #     print(f"[DEBUG] Running metadata command: {' '.join(meta_cmd)}")
    #     try:
    #         # subprocess.run(meta_cmd, check=True)
    #         with open(meta_path, "w") as meta_file:
    #             subprocess.run(meta_cmd, check=True, stdout=meta_file)
    #         print(f"Metadata file {sra} downloaded to {sra_dir}")
    #         obj["quantifications"][name]["sra"][sra]["meta_path"] = meta_path
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error while downloading Metadata file {sra}: {e}")
    #         return

    # print(f"SRA prefetch Done for {sra_list}")

    

    # ## Once downloaded, perform fastq-dump on the corresponding SRA files
    # for sra in sra_list:
    #     print(f"[DEBUG] Running fastq-dump for SRA accession: {sra}")
    #     sra_path = f"./data/{obj_name}/{name}/{sra}/{sra}.sra"
    #     fastqdump_path = f"./data/{obj_name}/{name}/{sra}/"
    #     print(f"[DEBUG] SRA path: {sra_path}")
    #     print(f"[DEBUG] Fastq output directory: {fastqdump_path}")

    #     if paired:
    #         fastqdump_cmd = [
    #             "fastq-dump",
    #             "--split-files",
    #             "-O",
    #             fastqdump_path,
    #             sra_path
    #         ]
    #         print(f"[DEBUG] Running paired-end fastq-dump command: {' '.join(fastqdump_cmd)}")
    #     else:
    #         fastqdump_cmd = [
    #             "fastq-dump",
    #             "-O",
    #             fastqdump_path,
    #             sra_path
    #         ]
    #         print(f"[DEBUG] Running single-end fastq-dump command: {' '.join(fastqdump_cmd)}")
        
    #     try:
    #         # Uncomment the next line to enable fastq-dump execution
    #         subprocess.run(fastqdump_cmd, check=True)
    #         print(f"Fastq files generated for {sra}")
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error while running fastq-dump for {sra}: {e}")
    #         return

    # print(f"Fastq dump Done for {sra_list}")

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
        print(f"[DEBUG] Gene counts DataFrame shape: {gene_counts_df.shape}")
        print(f"[DEBUG] Gene counts DataFrame head:\n{gene_counts_df.head()}")
        #columns names are 
        print(f"[DEBUG] Gene counts DataFrame columns:\n{gene_counts_df.columns}")
        # Transpose to make genes as columns
        gene_counts_df.T.to_csv(gene_count_csv)

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
    
    ##debug output
    print(f"Object Name: {obj_name}")
    
    config = None
    with open('config.json', 'r') as f:
        config = json.load(f)

    obj = config['objects'][obj_name]
    
    # "path_to_gene_count"
    path_to_gene_count = obj["quantifications"][quant_name]["path_to_gene_count"]

    df = pd.read_csv(path_to_gene_count)
    array_with_genes = df.values
    
    # Fetch the row where the first column matches the target string
    gene_row = array_with_genes[array_with_genes[:, 0] == target_gene]
    values = gene_row.iloc[0, 1:].values
    labels = df.columns[1:]

    plt.bar(labels, values)

def plot_gene_abundances(args):
    obj_name = args.obj
    quant_name = args.quantification_name
    target_gene = args.gene
    
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
        df = pd.read_csv(path_to_gene_count)
        if df.empty:
            raise ValueError(f"The gene counts file '{path_to_gene_count}' is empty.")
        
        # Ensure the target gene exists in the dataframe
        if target_gene not in df.iloc[:, 2].values:
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
        plt.ylabel("Expression Level")
        plt.title(f"Gene Abundance for {target_gene}")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.show()

    except IndexError as e:
        print(f"[ERROR] Index error occurred while processing gene row: {e}")
    except ValueError as e:
        print(f"[ERROR] Value error: {e}")
    except Exception as e:
        print(f"[ERROR] An unexpected error occurred during plotting: {e}")


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
