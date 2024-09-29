from pytximport import tximport
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

inference = DefaultInference(n_cpus=8)
t2g_full = pd.read_csv("/data1/yashjonjale/t2g.txt", sep="\t",
                       names=["transcript_id", "gene_id", "gene_name", "transcript_name", "chrom", "start", "end", "strand"])
t2g = t2g_full[["transcript_id", "gene_id"]]



# Load metadata
metadata = pd.read_csv("/home/schoudhary/github/Ehuxleyi/SRP182772/SRP182772.tsv", sep="\t")
metadata.set_index("run_accession", inplace=True)
metadata = metadata.sort_values("library_name")
metadata["library_name"] = pd.Categorical(metadata["library_name"], categories=metadata["library_name"].unique())

files = [f"/home/schoudhary/github/Ehuxleyi/SRP182772/snakemake_out/counts/{acc}/abundance.tsv" for acc in metadata.index]
#target_id       length  eff_length      est_counts      tpm
result = tximport(file_paths=files, data_type='kallisto', transcript_gene_map=t2g,
               id_column = "target_id",
               counts_column = "est_counts",
               length_column = "eff_length", abundance_column = "tpm")

result.X = result.X.round().astype(int)
#result = result[:, result.X.max(axis=0) > 10].copy()

metadata["type"] = metadata["library_name"].str.split("-", n=1).str[0]
metadata["type"] = pd.Categorical(metadata["type"], categories=metadata["type"].unique())
metadata["replicate"] = metadata["library_name"].str.split("-", n=1).str[1]

result.obs["type"] = metadata["type"]
result.obs["replicate"] = metadata["replicate"]

import pyranges as pr
gtf = pr.read_gtf("/data1/yashjonjale/igem_stuff/ehux_study/sra_data/working_folder/kb_ref_out/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf")
gtf_gene = gtf.df[gtf.df["Feature"] == "gene"][["gene_id", "gene_name"]].drop_duplicates()
gtf_gene["gene_name"].fillna(gtf_gene["gene_id"], inplace=True)
gtf_gene.set_index("gene_id", inplace=True)


common_genes = result.var.index.intersection(gtf_gene.index)

result = result[:, result.var.index.isin(common_genes)]
result.obs = metadata
#.loc[result.obs.index]

dds = DeseqDataSet(adata=result, design_factors="type",     refit_cooks=True,
                       inference=DefaultInference(n_cpus=8),
                       quiet=True,
                   )
dds.deseq2()


#vsd = vst(dds)
#vsd_data = vsd.transform()
dds.vst(use_design=True)

from sklearn.decomposition import PCA
df = dds.layers["vst_counts"]
pca_vsd =  PCA(n_components=2)
pca_vsd.fit(df)
pca_vsd = pca_vsd.transform(df)

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
sns.scatterplot(x=pca_vsd[:, 0], y=pca_vsd[:, 1], hue=metadata["type"])
plt.title("PCA (VST) by Type")

plt.subplot(1, 2, 2)
sns.scatterplot(x=pca_vsd[:, 0], y=pca_vsd[:, 1], hue=metadata["replicate"])
plt.title("PCA (VST) by Replicate")

plt.tight_layout()
plt.show()

plt.savefig("test.pdf")
