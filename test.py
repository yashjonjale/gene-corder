from pytximport.utils import create_transcript_to_gene_map_from_gtf_annotation

transcript_gene_map_from_gtf = create_transcript_to_gene_map_from_gtf_annotation(
    "",
    target_field="gene_name",
)
# transcript_gene_map_from_gtf.head(5)

import numpy as np
import pandas as pd
from pytximport import tximport


txi = tximport(
    ["../../test/data/salmon/multiple/Sample_1.sf", "../../test/data/salmon/multiple/Sample_2.sf"],
    "salmon",
    transcript_gene_map_mouse,
    counts_from_abundance="length_scaled_tpm",
    output_type="xarray",  # or "anndata"
)
# txi

from pytximport.utils import replace_transcript_ids_with_names

transcript_name_map_human = create_transcript_to_gene_map("human", target_field="external_transcript_name")
# transcript_name_map_human.head(5)


from pytximport.utils import replace_gene_ids_with_names




