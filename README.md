# GENECorder

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Getting Started](#getting-started)
  - [Basic Commands](#basic-commands)
- [Advanced Features](#advanced-features)
  - [Integration with Custom Databases](#integration-with-custom-databases)
    - [Example Workflow with Custom Genome](#example-workflow-with-custom-genome)
- [Integrated Software Tools](#integrated-software-tools)
- [Support and Contributions](#support-and-contributions)
- [License](#license)

---

## Overview

The **RNASeq Data Analyser** is an open-source, easy-to-use tool designed to help researchers analyze RNA sequencing (RNA-seq) datasets, regardless of their programming or bioinformatics expertise. Built to support bulk RNA-seq data from organisms listed on Ensembl, the tool simplifies data processing, quantifies gene expression, and produces an easily interpretable output file (`gene_counts.csv`) for downstream analysis. It integrates multiple bioinformatics tools behind the scenes, offering a smooth and automated workflow.

---

## Key Features

- **Easy Integration with Public Data**: Automatically fetches genomic and transcriptomic data from the Ensembl database, supporting a wide range of organisms.
- **Built on Proven Bioinformatics Tools**: Incorporates trusted bioinformatics software like `kallisto`, `kb-python`, and `tximport` for high-quality results.
- **User-Friendly Command-Line Interface**: Simple commands provided for easy interaction without requiring programming knowledge.
- **Cross-Database Compatibility**: Can map gene indices from multiple databases, enabling analysis with custom genomes and annotations.
- **Visual and Statistical Insights**: Built-in functions for gene expression plots and co-expression analysis.

---

## System Requirements

- **Supported Platforms**: macOS and Linux.
- **Required Dependencies**:
  - **Anaconda**: For managing Python dependencies.
  - **R**: Required for `tximport` and `DESeq2` integration.
  - **Internet Access**: Needed to download required genome and transcriptome data from Ensembl.
  - **Disk Space**: Ensure sufficient space is available for storing downloaded genome files.

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

## Getting Started

### Basic Commands

1. **List available organisms**:

   ```bash
   ./main.py list
   ```

2. **Instantiate an organism**:

   ```bash
   ./main.py instantiate --organism "Homo_sapiens"
   ```

3. **Quantize RNA-Seq data**:

   ```bash
   ./main.py quantize --sra "SRR123456" --name "my_quant"
   ```

4. **List quantifications**:

   ```bash
   ./main.py list_quant
   ```

5. **Plot gene abundances**:

   ```bash
   ./main.py plot_gene_abundances --gene "BRCA1"
   ```

6. **Generate correlation matrix**:

   ```bash
   ./main.py generate_correlation_matrix --genes "BRCA1,TP53"
   ```

7. **Perform DESeq2 analysis**:

   ```bash
   ./main.py deseq_analyse --gene "BRCA1"
   ```

---

## Advanced Features

### Integration with Custom Databases

The **RNASeq Data Analyser** allows you to integrate custom genomes and annotations that are not available in the Ensembl database. This is particularly useful when working with specialized organisms, strains, or when you have custom genome assemblies. The tool uses **Multi-Liftoff** under the hood to map annotations from Ensembl to your custom genome, enabling seamless analysis with your data.

#### Command: `intersect`

**Description**:

The `intersect` command creates a mapping between Ensembl data and your custom genome and annotations.

**Usage**:

```bash
./main.py intersect --organism "Organism_Name" --custom_genome "path/to/custom_genome.fa" --custom_gtf "path/to/custom_annotations.gtf"
```

- `--organism`: The name of the organism you've instantiated.
- `--custom_genome`: Path to your custom genome FASTA file.
- `--custom_gtf`: Path to your custom GTF annotation file.

**What Happens Under the Hood**:

- **Annotation Mapping**: Uses Multi-Liftoff to map Ensembl annotations to your custom genome.
- **Mapping File Generation**: Creates a mapping file linking Ensembl gene/transcript IDs to your custom IDs.
- **Configuration Update**: Updates `config.json` to include the custom mappings.
- **Downstream Integration**: Allows you to use gene indices from your custom database in all downstream analyses.

---

#### Example Workflow with Custom Genome

**Scenario**:

You are working with a custom strain of *Saccharomyces cerevisiae* that has a slightly different genome assembly from the one available in Ensembl. You want to analyze RNA-Seq data using your custom genome and annotations.

**Steps**:

1. **Instantiate the Organism**:

   Even though your custom strain isn't in Ensembl, instantiate the closest available organism.

   ```bash
   ./main.py instantiate --organism "Saccharomyces_cerevisiae"
   ```

2. **Run the Intersect Command**:

   Map Ensembl annotations to your custom genome.

   ```bash
   ./main.py intersect --organism "Saccharomyces_cerevisiae" --custom_genome "/path/to/custom_genome.fa" --custom_gtf "/path/to/custom_annotations.gtf"
   ```

   *Explanation*:

   - **Custom Genome and GTF**: Provide paths to your custom genome and GTF files.
   - **Under the Hood**: The tool runs Multi-Liftoff to map Ensembl annotations to your custom genome, generating a mapping file.

3. **Quantize RNA-Seq Data with Custom Genome**:

   Fetch and quantify RNA-Seq data using your custom genome.

   ```bash
   ./main.py quantize --sra "SRR654321" --name "custom_quant"
   ```

   *Note*: The tool uses the custom genome and mappings for quantification.

4. **Plot Gene Expression with Custom Indices**:

   You can now query genes using indices from your custom annotations.

   ```bash
   ./main.py plot_gene_abundances --gene "YFG1_custom"
   ```

   *Explanation*:

   - **Custom Gene ID**: Replace `YFG1_custom` with your gene ID from the custom annotations.
   - **Under the Hood**: The tool translates the custom gene ID to the corresponding Ensembl ID using the mapping file.

5. **Generate Correlation Matrix with Custom Genes**:

   ```bash
   ./main.py generate_correlation_matrix --genes "YFG1_custom,YFG2_custom"
   ```

   *Explanation*:

   - **Custom Gene IDs**: Provide a comma-separated list of custom gene IDs.
   - **Analysis**: The tool performs co-expression analysis using your custom genes.

---

#### Additional Examples

**Example 1: Working with a Plant Genome Not in Ensembl**

Let's say you're studying a plant species, *Plantus exampleus*, not available in Ensembl.

1. **Instantiate a Close Relative**:

   ```bash
   ./main.py instantiate --organism "Arabidopsis_thaliana"
   ```

2. **Intersect with Custom Genome**:

   ```bash
   ./main.py intersect --organism "Arabidopsis_thaliana" --custom_genome "/path/to/plantus_exampleus_genome.fa" --custom_gtf "/path/to/plantus_exampleus_annotations.gtf"
   ```

3. **Proceed with Analysis**:

   - Quantize RNA-Seq data.
   - Use custom gene IDs in your analyses.

**Example 2: Integrating a Disease-Specific Database**

Suppose you have specialized annotations for a disease model in mice.

1. **Instantiate Mouse Genome**:

   ```bash
   ./main.py instantiate --organism "Mus_musculus"
   ```

2. **Intersect with Custom Annotations**:

   ```bash
   ./main.py intersect --organism "Mus_musculus" --custom_genome "/path/to/custom_mouse_genome.fa" --custom_gtf "/path/to/disease_specific_annotations.gtf"
   ```

3. **Analyze Using Disease-Specific Indices**:

   - Use gene indices from the disease-specific database.
   - Perform quantification, plotting, and differential expression analysis.

---

## Integrated Software Tools

| Tool            | Functionality                                        | Purpose in RNASeq Data Analyser                   |
|-----------------|------------------------------------------------------|---------------------------------------------------|
| **Kallisto**    | RNA-seq transcript quantification                    | Processes RNA-seq reads into transcript abundances |
| **kb-python**   | RNA-seq preprocessing, including handling FASTQ files | Fetches, prepares, and processes RNA-seq data     |
| **Minimap2**    | Sequence alignment                                   | Aligns sequences for mapping genomic features     |
| **Multi-Liftoff** | Cross-database mapping and alignment                 | Maps gene indices across multiple databases       |
| **Samtools**    | Manipulation of alignment data (SAM/BAM files)        | Manages sequence alignments and data formats      |
| **Bamtools**    | Utilities for BAM file manipulation                   | Manages BAM files, used for storing aligned sequences |
| **Bedtools**    | Operates on genomic intervals (BED files)             | Processes and intersects data with genome annotations |
| **Gffread**     | GTF/GFF3 file processing                             | Extracts and converts gene annotations            |
| **Biopython**   | Biological computation and sequence manipulation     | Provides utilities for biological data processing |
| **EMBOSS**      | Bioinformatics tools for sequence analysis           | Adds support for additional sequence operations   |
| **Ensembl API** | Interface for accessing genomic data from Ensembl    | Downloads necessary organism-specific data        |
| **Pydeseq2**    | Differential expression analysis using DESeq2        | Calculates gene expression differences            |
| **Tximport**    | Converts transcript abundances to gene-level counts  | Handles conversion of quantified transcripts      |
| **R and DESeq2**| Statistical computing and differential expression analysis | Performs normalization and statistical testing |

---

## Support and Contributions

If you have any questions or encounter issues, feel free to reach out to the development team. We welcome contributions—whether through suggestions, bug reports, or code improvements. You can contribute by submitting pull requests on GitHub or by contacting the team directly.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

# Important Note

The code examples and commands provided assume that all necessary dependencies and environment configurations are correctly set up. Ensure that you have the required tools installed and accessible in your system's PATH.
