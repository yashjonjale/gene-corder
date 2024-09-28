import argparse
import os
import sys
import json
import subprocess
import requests
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def save_config(config_json):
    ##save the config file
    with open('config.json', 'w') as f:
        json.dump(config_json, f, indent=4)


def save_config(config):
    with open('config.json', 'w') as f:
        json.dump(config, f, indent=4)



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

    # Directory to store outputs
    output_dir = f"./data/{name}"
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Index genome using kallisto index
    kallisto_index_path = os.path.join(output_dir, f"{name}_kallisto.idx")
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

    # # Step 2: Generate t2g.tsv using kb ref
    # t2g_path = os.path.join(output_dir, f"{name}_t2g.tsv")
    # print(f"Generating t2g.tsv using kb ref...")
    # kb_ref_cmd = [
    #     "kb", "ref", "-i", kallisto_index_path,
    #     "-g", gtf_path,
    #     "-f1", transcript_path,
    #     "-o", output_dir
    # ]
    
    # try:
    #     subprocess.run(kb_ref_cmd, check=True)
    #     print(f"t2g.tsv generated at {t2g_path}")
    #     obj["index_paths"]["t2g"] = t2g_path
    # except subprocess.CalledProcessError as e:
    #     print(f"Error while running kb ref: {e}")
    #     return

    # # Step 3: Save the updated config file with the new paths
    # config['objects'][name] = obj
    # save_config(config)
    # print(f"Updated config file with new paths for {name}.")
    
    # print(f"Genome indexing and t2g generation complete. Ready for downstream quantification!")

