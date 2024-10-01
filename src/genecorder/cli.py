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
    print("DSADASDA")
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
    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    if args.command == "instantiate":
        instantiate(args)
    elif args.command == "quantize":
        quantize(args)



if __name__ == "__main__":
    parse_args(sys.argv[1:])
