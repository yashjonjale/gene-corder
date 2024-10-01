"""Console script for genecorder."""

import argparse
import os
import re
import sys
import warnings
from io import StringIO
from textwrap import dedent

import pandas as pd

from . import __version__
from .genecorder import *


class ArgParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


class CustomFormatterArgP(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def create_config():
    if not os.path.exists("config.json"):
        print("[DEBUG] config.json does not exist. Creating a new one.")
        with open("config.json", "w") as f:
            json.dump({"objects": {}}, f, indent=4)
    else:
        print("[DEBUG] config.json found.")


def parse_args(args=None):
    parser = ArgParser(
        description="genecorder (v{})".format(__version__),
        formatter_class=CustomFormatterArgP,
    )
    subparsers = parser.add_subparsers(title="subcommands", dest="command")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=__version__),
    )

    # Instantiate
    parser_instantiate = subparsers.add_parser(
        "instantiate",
        help="Instantiate an object of an organism, providing paths to the genome, transcriptome, and annotation files, ",
    )
    parser_instantiate.add_argument("--organism", required=True, help="Organism name")
    parser_instantiate.add_argument(
        "--name", required=True, help="Name for this quantification instance"
    )
    parser_instantiate.add_argument(
        "--desc", required=False, help="Description for this quantification"
    )
    parser_instantiate.add_argument(
        "--transcriptome_path",
        required=False,
        help="Path to transcriptome file in fasta format",
    )
    parser_instantiate.add_argument(
        "--genome_path", required=True, help="Path to genome file in fasta format"
    )
    parser_instantiate.add_argument(
        "--gtf_path",
        required=True,
        help="Path to the annotation file in GTF/GFF format",
    )
    parser_instantiate.set_defaults(func=instantiate_organism)

    # Quantize
    parser_quantize = subparsers.add_parser("quantize", help="Quantize RNA-Seq data")
    parser_quantize.add_argument(
        "--sra", required=True, help="Comma-separated SRA accession codes"
    )
    parser_quantize.add_argument(
        "--name", required=True, help="Name for this quantification"
    )
    parser_quantize.add_argument("--obj", required=True, help="Object name")
    parser_quantize.add_argument(
        "--paired", required=False, action="store_true", help="Paired-end reads"
    )
    parser_quantize.set_defaults(func=quantize)

    # Uncomment and implement additional subparsers as needed
    parser_list_quant = subparsers.add_parser("listquants", help="List quantifications")
    parser_list_quant.add_argument("--obj", required=True, help="Object name")
    parser_list_quant.set_defaults(func=list_quantifications)

    parser_plot = subparsers.add_parser("plotga", help="Plot gene abundances")
    parser_plot.add_argument("--gene", required=True, help="Gene name")
    parser_plot.add_argument(
        "--named", required=False, action="store_true", help="Gene name is provided"
    )
    parser_plot.add_argument("--obj", required=True, help="Object name")
    parser_plot.add_argument(
        "--quantification_name", required=True, help="Quantification name"
    )
    parser_plot.add_argument("--output", required=True, help="Output file path")
    parser_plot.set_defaults(func=plot_gene_abundances)

    parser_corr = subparsers.add_parser("corr", help="Generate correlation matrix")
    parser_corr.add_argument(
        "--genes", required=True, help="newline separated list of genes"
    )
    parser_corr.add_argument("--obj", required=True, help="Object name")
    parser_corr.add_argument(
        "--quantification_name", required=True, help="Quantification name"
    )
    parser_corr.add_argument("--output_dir", required=True, help="Output directory")
    parser_corr.set_defaults(func=generate_correlation_matrix)

    parser_name2id = subparsers.add_parser(
        "name2id", help="Convert gene name to gene ID"
    )
    parser_name2id.add_argument("--gene_name", required=True, help="Gene name")
    parser_name2id.add_argument("--obj", required=True, help="Object name")

    parser_name2id.set_defaults(func=name2id)

    parser_gene2fasta = subparsers.add_parser(
        "gene2fasta", help="Extract gene sequence to FASTA"
    )
    parser_gene2fasta.add_argument("--gene", required=True, help="Gene name or ID")
    parser_gene2fasta.add_argument("--obj", required=True, help="Object name")
    parser_gene2fasta.add_argument(
        "--output_dir", required=True, help="Output directory"
    )
    parser_gene2fasta.add_argument(
        "--named", required=False, action="store_true", help="Gene name is provided"
    )

    parser_gene2fasta.set_defaults(func=gene2fasta)

    parser_list_objs = subparsers.add_parser("listobjs", help="List objects")
    parser_list_objs.set_defaults(func=list_objs)

    parser_deseq = subparsers.add_parser("deseq", help="Perform DESeq2 analysis")
    parser_deseq.add_argument("--obj", required=True, help="Object name")
    parser_deseq.add_argument(
        "--quantification_name", required=True, help="Quantification name"
    )
    parser_deseq.add_argument(
        "--srp", required=True, help="Comma-separated SRA accession codes"
    )
    parser_deseq.add_argument(
        "--paired", required=False, action="store_true", help="Paired-end reads"
    )
    parser_deseq.add_argument("--output_dir", required=True, help="Output directory")

    parser_deseq.set_defaults(func=quant_deseq)

    parser_remove = subparsers.add_parser("remove", help="Remove object")
    parser_remove.add_argument("--obj", required=True, help="Object name")
    parser_remove.set_defaults(func=remove_object)

    args = parser.parse_args()
    print(f"[DEBUG] Parsed arguments: {args}")

    if hasattr(args, "func"):
        print(f"[DEBUG] Executing command: {args.command}")
        create_config()
        args.func(args)
    else:
        print("[DEBUG] No command provided. Displaying help.")
        parser.print_help()


if __name__ == "__main__":
    parse_args(sys.argv[1:])
