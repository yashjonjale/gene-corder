# GENECorder

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Integrated Software Tools](#integrated-software-tools)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Commands and Getting Started](#commands-and-getting-started)
- [Examples](#example-workflows)
- [Additional Note](#additional-notes)
- [Support and Contributions](#support-and-contributions)
- [License](#license)

---


## Overview

The **GENECorder** is an open-source, easy-to-use tool designed to help researchers analyze RNA sequencing (RNA-seq) datasets, regardless of their programming or bioinformatics expertise. Built to support bulk RNA-seq data from organisms listed on Ensembl, the tool simplifies data processing, quantifies gene expression, and produces an easily interpretable output file (`gene_counts.csv`) for downstream analysis. It integrates multiple bioinformatics tools behind the scenes, offering a smooth and automated workflow.


---

## Key Features

- **Easy Integration with Public Data**: Automatically fetches genomic and transcriptomic data from the Ensembl database, supporting a wide range of organisms.
- **Built on Proven Bioinformatics Tools**: Incorporates trusted bioinformatics software like `kallisto`, `kb-python`, and `pytximport` for high-quality results.
- **User-Friendly Command-Line Interface**: Simple commands provided for easy interaction without requiring programming knowledge.
- **Cross-Database Compatibility**: Can map gene indices from multiple databases, enabling analysis with custom genomes and annotations.
- **Visual and Statistical Insights**: Built-in functions for gene expression plots and co-expression analysis.

---

## Integrated Software Tools

| Tool            | Functionality                                        | Purpose in RNASeq Data Analyser                   |
|-----------------|------------------------------------------------------|---------------------------------------------------|
| **Kallisto**    | RNA-seq transcript quantification                    | Processes RNA-seq reads into transcript abundances |
| **kb-python**   | RNA-seq preprocessing, including handling FASTQ files | Fetches, prepares, and processes RNA-seq data     |
| **Samtools**    | Manipulation of alignment data (SAM/BAM files)        | Manages sequence alignments and data formats      |
| **Bamtools**    | Utilities for BAM file manipulation                   | Manages BAM files, used for storing aligned sequences |
| **Bedtools**    | Operates on genomic intervals (BED files)             | Processes and intersects data with genome annotations |
| **Gffread**     | GTF/GFF3 file processing                             | Extracts and converts gene annotations            |
| **Biopython**   | Biological computation and sequence manipulation     | Provides utilities for biological data processing |
| **Pydeseq2**    | Differential expression analysis using DESeq2        | Calculates gene expression differences            |
| **PyTximport**    | Converts transcript abundances to gene-level counts  | Handles conversion of quantified transcripts      |


---

## System Requirements

- **Supported Platforms**: macOS and Linux.
- **Required Dependencies**:
  - **Anaconda**: For managing Python dependencies. Nothing else needed, we take care of the rest.
  - **Internet Access**: Needed to prefetch SRA files from NCBI
  - **Disk Space**: Ensure sufficient space is available for storing downloaded genome files which depends on the extent of your analyses. Guven it is bioinformatics throughput, memory requirement may sometimes be over 100 GBs

---

## Installation

1. **Clone the GitHub repository**:

   ```bash
   git clone <repository-url>
   ```

2. **Ensure Anaconda is installed**:

   Download it from [here](https://www.anaconda.com/products/individual).

3. **Run the `make` command**:

   ```bash
   make
   ```

   This command will install all the necessary dependencies and set up the environment automatically.

---


## Usage

The general usage pattern for GENECorder is:

```bash
python main.py <command> [options]
```

To see the help message and list of available commands:

```bash
python main.py --help
```


---

## Commands and Getting Started

### 1. `instantiate`

**Description:** Instantiate an object for an organism. This sets up the necessary indices and mappings for downstream analyses.

**Usage:**

```bash
python main.py instantiate --organism <organism_name> --name <object_name> [--desc <description>] [--transcriptome_path <transcriptome.fasta>] --genome_path <genome.fasta> --gtf_path <annotation.gtf>
```

**Arguments:**

- `--organism`: **(Required)** Name of the organism (e.g., `Homo sapiens`).
- `--name`: **(Required)** A unique name for this object (used in future commands).
- `--desc`: **(Optional)** Description for this object.
- `--transcriptome_path`: **(Optional)** Path to transcriptome FASTA file. If not provided, it will be generated from the genome and GTF files.
- `--genome_path`: **(Required)** Path to genome FASTA file.
- `--gtf_path`: **(Required)** Path to annotation GTF/GFF file.

**Example:**

```bash
python main.py instantiate --organism "Homo sapiens" --name human_obj --genome_path /path/to/human_genome.fa --gtf_path /path/to/human_annotation.gtf
```

---

### 2. `quantize`

**Description:** Quantify RNA-Seq data using Kallisto for the specified object.

**Usage:**

```bash
python main.py quantize --sra <sra_accession_codes> --name <quantification_name> --obj <object_name> [--paired]
```

**Arguments:**

- `--sra`: **(Required)** Comma-separated list of SRA accession codes (e.g., `SRR1234567,SRR1234568`).
- `--name`: **(Required)** Name for this quantification.
- `--obj`: **(Required)** Name of the object created with `instantiate`.
- `--paired`: **(Optional)** Include if reads are paired-end.

**Example:**

```bash
python main.py quantize --sra SRR1234567,SRR1234568 --name quant1 --obj human_obj --paired
```

---

### 3. `list_quant`

**Description:** List all quantifications associated with a specified object.

**Usage:**

```bash
python main.py list_quant --obj <object_name>
```

**Arguments:**

- `--obj`: **(Required)** Name of the object.

**Example:**

```bash
python main.py list_quant --obj human_obj
```

---

### 4. `plot_gene_abundances`

**Description:** Plot the abundances of a specified gene across samples in a quantification.

**Usage:**

```bash
python main.py plot_gene_abundances --gene <gene_name_or_id> [--named] --obj <object_name> --quantification_name <quantification_name> --output <output_file>
```

**Arguments:**

- `--gene`: **(Required)** Gene name or ID to plot.
- `--named`: **(Optional)** Include if providing a gene name instead of an ID.
- `--obj`: **(Required)** Name of the object.
- `--quantification_name`: **(Required)** Name of the quantification.
- `--output`: **(Required)** Path to save the output plot (e.g., `abundance_plot.png`).

**Example:**

Plotting by gene name:

```bash
python main.py plot_gene_abundances --gene BRCA1 --named --obj human_obj --quantification_name quant1 --output brca1_abundance.png
```

Plotting by gene ID:

```bash
python main.py plot_gene_abundances --gene ENSG00000012048 --obj human_obj --quantification_name quant1 --output brca1_abundance.png
```

---

### 5. `generate_correlation_matrix`

**Description:** Generate a correlation matrix and heatmap for a list of genes based on their expression levels.

**Usage:**

```bash
python main.py generate_correlation_matrix --genes <gene_list_file> --obj <object_name> --quantification_name <quantification_name> --output_dir <output_directory>
```

**Arguments:**

- `--genes`: **(Required)** Path to a text file containing a list of gene IDs or names (one per line).
- `--obj`: **(Required)** Name of the object.
- `--quantification_name`: **(Required)** Name of the quantification.
- `--output_dir`: **(Required)** Directory to save the correlation matrix and heatmap.

**Example:**

```bash
python main.py generate_correlation_matrix --genes genes.txt --obj human_obj --quantification_name quant1 --output_dir ./correlation_results/
```

---

### 6. `name2id`

**Description:** Convert a gene name to its corresponding gene ID.

**Usage:**

```bash
python main.py name2id --gene_name <gene_name> --obj <object_name>
```

**Arguments:**

- `--gene_name`: **(Required)** The gene name to convert.
- `--obj`: **(Required)** Name of the object.

**Example:**

```bash
python main.py name2id --gene_name BRCA1 --obj human_obj
```

---

### 7. `gene2fasta`

**Description:** Extract the sequence of a gene and save it as a FASTA file.

**Usage:**

```bash
python main.py gene2fasta --gene <gene_name_or_id> --obj <object_name> --output_dir <output_directory> [--named]
```

**Arguments:**

- `--gene`: **(Required)** Gene name or ID to extract.
- `--obj`: **(Required)** Name of the object.
- `--output_dir`: **(Required)** Directory to save the FASTA file.
- `--named`: **(Optional)** Include if providing a gene name instead of an ID.

**Example:**

Extracting by gene name:

```bash
python main.py gene2fasta --gene BRCA1 --named --obj human_obj --output_dir ./gene_sequences/
```

Extracting by gene ID:

```bash
python main.py gene2fasta --gene ENSG00000012048 --obj human_obj --output_dir ./gene_sequences/
```

---

### 8. `list_objs`

**Description:** List all instantiated objects.

**Usage:**

```bash
python main.py list_objs
```

**Example:**

```bash
python main.py list_objs
```

---

### 9. `deseq_analyse`

**Description:** Perform differential expression analysis using DESeq2 on RNA-Seq data from a specified SRA project.

**Usage:**

```bash
python main.py deseq_analyse --obj <object_name> --quantification_name <analysis_name> --srp <sra_project_code> [--paired] --output_dir <output_directory>
```

**Arguments:**

- `--obj`: **(Required)** Name of the object.
- `--quantification_name`: **(Required)** Name for this analysis.
- `--srp`: **(Required)** SRA project accession code (e.g., `SRP123456`).
- `--paired`: **(Optional)** Include if reads are paired-end.
- `--output_dir`: **(Required)** Directory to save analysis results.

**Example:**

```bash
python main.py deseq_analyse --obj human_obj --quantification_name de_analysis1 --srp SRP123456 --paired --output_dir ./deseq_results/
```

---

### 10. `remove`

**Description:** Remove an instantiated object and its associated data.

**Usage:**

```bash
python main.py remove --obj <object_name>
```

**Arguments:**

- `--obj`: **(Required)** Name of the object to remove.

**Example:**

```bash
python main.py remove --obj human_obj
```

---

## Additional Notes

- Ensure that you have a stable internet connection when running commands that download data (e.g., `quantize`, `deseq_analyse`).
- All data and outputs are organized under the `data/` directory, structured by object and quantification names.
- The tool maintains a `config.json` file to track objects and their associated data. Do not modify this file manually.
- Use the `--paired` flag if your RNA-Seq data is from paired-end sequencing; otherwise, omit it for single-end data.
- For the `generate_correlation_matrix` command, the gene list file should contain one gene name or ID per line.

---

