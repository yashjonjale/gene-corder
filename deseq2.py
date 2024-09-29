from pytximport import tximport
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

inference = DefaultInference(n_cpus=8)
t2g_full = pd.read_csv("/data1/yashjonjale/igem_stuff/ehux_study/sra_data/working_folder/w_datasets2/t2g.txt", sep="\t",
                       names=["txid", "gene_id", "gene_name", "chrom", "start", "end", "strand"])
t2g = t2g_full[["txid", "gene_id"]]

# Load metadata
metadata = pd.read_csv("/home/schoudhary/github/Ehuxleyi/SRP182772/SRP182772.tsv", sep="\t")
metadata.set_index("run_accession", inplace=True)
metadata = metadata.sort_values("library_name")
metadata["library_name"] = pd.Categorical(metadata["library_name"], categories=metadata["library_name"].unique())

files = [f"/home/schoudhary/github/Ehuxleyi/SRP182772/snakemake_out/counts/{acc}/abundance.tsv" for acc in metadata.index]
txi = pytximport.import_kallisto(files=files, tx2gene=t2g)
counts = txi["counts"].round(0)

import pyranges as pr
gtf = pr.read_gtf("/data1/yashjonjale/igem_stuff/ehux_study/sra_data/working_folder/kb_ref_out/Emiliania_huxleyi.Emiliana_huxleyi_CCMP1516_main_genome_assembly_v1.0.59.gtf")
gtf_gene = gtf.df[gtf.df["Feature"] == "gene"][["gene_id", "gene_name"]].drop_duplicates()
gtf_gene["gene_name"].fillna(gtf_gene["gene_id"], inplace=True)
gtf_gene.set_index("gene_id", inplace=True)

common_genes = counts.index.intersection(gtf_gene.index)
counts = counts.loc[common_genes]
counts.index = pd.Index(pd.Series(common_genes).apply(lambda x: pd.unique(gtf_gene.loc[x, "gene_name"])))

metadata = metadata.loc[counts.columns]
metadata["type"] = metadata["library_name"].str.split("-", n=1).str[0]
metadata["type"] = pd.Categorical(metadata["type"], categories=metadata["type"].unique())
metadata["replicate"] = metadata["library_name"].str.split("-", n=1).str[1]

dds = DeseqDataSet(counts=counts, clinical=metadata, design_factors="type")
dds.fit()

vsd = VST(dds)
vsd_data = vsd.transform()

rld = RLOG(dds)
rld_data = rld.transform()

pca_vsd = vsd.plot_pca(group="type")
pca_replicate = vsd.plot_pca(group="replicate")

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
sns.scatterplot(x=pca_vsd[:, 0], y=pca_vsd[:, 1], hue=metadata["type"])
plt.title("PCA (VST) by Type")

plt.subplot(1, 2, 2)
sns.scatterplot(x=pca_replicate[:, 0], y=pca_replicate[:, 1], hue=metadata["replicate"])
plt.title("PCA (VST) by Replicate")

plt.tight_layout()
plt.show()


