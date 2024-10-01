# GENECorder CLI Tool Documentation

## Introduction

**GENECorder** is a command-line interface (CLI) tool designed to streamline and automate various genomic data analysis tasks, particularly focusing on RNA-Seq data quantification, gene expression analysis, and sequence extraction. It provides a user-friendly interface to perform common bioinformatics workflows, integrating tools like Kallisto, DESeq2, and others to facilitate gene quantification, differential expression analysis, and more.

With GENECorder, you can:

- Instantiate organism-specific objects to manage data and analyses.
- Quantify RNA-Seq data using Kallisto.
- Perform differential expression analysis with DESeq2.
- Plot gene abundances across samples.
- Generate correlation matrices for gene expression.
- Extract gene sequences into FASTA files.
- Convert between gene names and IDs.
- Manage and list analysis objects.

---

## Installation

To use GENECorder, ensure you have Python 3 installed along with the necessary dependencies. The required Python packages include:

- `argparse`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `pyranges`
- `pyfaidx`
- `Bio`
- `pytximport`
- `pydeseq2`
- `scikit-learn`

You may also need to install external bioinformatics tools that GENECorder relies on, such as:

- **Kallisto**: For RNA-Seq quantification.
- **SRA Toolkit**: For downloading sequencing data from NCBI SRA.
- **GFFRead**: For converting GTF files to FASTA.

Install the required Python packages using `pip`:

```bash
pip install pandas numpy matplotlib seaborn pyranges pyfaidx biopython pytximport pydeseq2 scikit-learn
```

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

## Commands

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

## Example Workflow

Below is an example of a typical workflow using GENECorder:

1. **Instantiate an Object for an Organism:**

   ```bash
   python main.py instantiate --organism "Homo sapiens" --name human_obj --genome_path /path/to/human_genome.fa --gtf_path /path/to/human_annotation.gtf
   ```

2. **Quantify RNA-Seq Data:**

   ```bash
   python main.py quantize --sra SRR1234567,SRR1234568 --name quant1 --obj human_obj --paired
   ```

3. **List Quantifications:**

   ```bash
   python main.py list_quant --obj human_obj
   ```

4. **Plot Gene Abundances:**

   ```bash
   python main.py plot_gene_abundances --gene BRCA1 --named --obj human_obj --quantification_name quant1 --output brca1_abundance.png
   ```

5. **Generate Correlation Matrix:**

   ```bash
   python main.py generate_correlation_matrix --genes genes.txt --obj human_obj --quantification_name quant1 --output_dir ./correlation_results/
   ```

6. **Perform Differential Expression Analysis:**

   ```bash
   python main.py deseq_analyse --obj human_obj --quantification_name de_analysis1 --srp SRP123456 --paired --output_dir ./deseq_results/
   ```

7. **Extract Gene Sequence:**

   ```bash
   python main.py gene2fasta --gene BRCA1 --named --obj human_obj --output_dir ./gene_sequences/
   ```

8. **Remove an Object:**

   ```bash
   python main.py remove --obj human_obj
   ```

---

## Troubleshooting

- **Config File Not Found:** Ensure you have run the `instantiate` command to create the initial configuration.
- **File Not Found Errors:** Verify that all file paths provided are correct and accessible.
- **Gene Not Found:** Ensure gene names or IDs are correct and present in the organism's annotation.
- **Internet Connection:** Some commands require downloading data from NCBI SRA; ensure you have internet access.

---

## Contact and Support

For issues, questions, or contributions, please refer to the project's GitHub repository or contact the maintainer.

---

*This documentation provides an overview of GENECorder's capabilities and usage. For detailed information, please consult the source code or reach out to the development team.*