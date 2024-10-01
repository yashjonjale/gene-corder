import gffutils
import sys
import subprocess
import os
import pyranges as pr
from collections import defaultdict
import gzip
import pandas as pd
import re


GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

def gtf_to_dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def write_gtf(df, output_file):
    with open(output_file, 'w') as f:
        for index, row in df.iterrows():
            attributes = []
            attributes.append(f'gene_name "{row["gene_name"]}";')
            if pd.notnull(row['transcript_id']):
                attributes.append(f'transcript_id "{row["transcript_id"]}";')
            if pd.notnull(row['protein_id']):
                attributes.append(f'protein_id "{row["protein_id"]}";')
            if pd.notnull(row['exon_number']):
                attributes.append(f'exon_number "{row["exon_number"]}";')
            if pd.notnull(row['frame']):
                attributes.append(f'frame "{row["frame"]}";')

            attributes_str = ' '.join(attributes)

            # Write the GTF line
            gtf_line = (
                f'{row["seqname"]}\t'
                f'{row["source"]}\t'
                f'{row["feature"]}\t'
                f'{int(row["start"])}\t'
                f'{int(row["end"])}\t'
                f'{row["score"] if pd.notnull(row["score"]) else "."}\t'
                f'{row["strand"]}\t'
                f'{row["frame"] if pd.notnull(row["frame"]) else "."}\t'
                f'{attributes_str}'
            )
            f.write(gtf_line + '\n')

def do_liftoff(query_fasta, target_fasta, input_gtf, output_gtf, format_input_gtf=False, output_gtf_db=None):
    if not output_gtf_db:
        output_gtf_db = "{}_gffutils.db".format(input_gtf)
    #f formatinput_gtf:

    featuredb = gffutils.create_db(
        data=input_gtf,
        dbfn=output_gtf_db,
        merge_strategy="create_unique",
        force=True,
        disable_infer_transcripts=True,
        disable_infer_genes=True,
        verbose=True,
    )
    flank_length = "0.1"

    overlap_s = "0.5"
    overlap_a = "0.5"

    threads = "4"
    cmd = ["liftoff",  "-db",  output_gtf_db, "-flank",  flank_length, "-s", overlap_s, "-a", overlap_a, "-o", output_gtf, "-p{}".format(threads), query_fasta, target_fasta]
    cmd = " ".join(cmd)
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

def compare_gtf(query_gtf, target_gtf):
    query = pr.read_gtf(query_gtf)
    target = pr.read_gtf(target_gtf)

    query_gene = query[query.Feature == "gene"]
    target_gene = target[target.Feature == "gene"]

    intersection = query_gene.intersect(target_gene)

    print(intersection)






if __name__ == "__main__":
    do_liftoff("/home/schoudhary/github/Ehuxleyi/genome/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.dna.toplevel.fa", 
               "/data1/yashjonjale/from_yash_local/Emihu1/Emihu1_AssemblyScaffolds_2021-09-22.fasta", 
               #"/home/schoudhary/github/Ehuxleyi/genome/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf",
               #"/data1/yashjonjale/from_yash_local/Emihu1/Emihu1_all_genes.gff", 
               "Emihu1/test.gtf", 
               "outputx_gtf.db")
    # compare_gtf("output.gtf", "/data1/yashjonjale/run_liftoff_arab/Arabidopsis_thaliana.TAIR10.59.gtf")
