import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from core.utils import get_sig_lfcs


SERVER_ROOT = Path("/Volumes/kml29j/snRNA_bench")
# PROJECT_ROOT = Path.cwd().parent
data_dir = SERVER_ROOT / "data"

h5ad_dir = data_dir / 'h5ad_files' / 'processed_laptop_sept26'


# Collect gene symbol → Ensembl ID mappings from all samples
# # sb_samples = {
#     "AF007_SB_SC1_A": "UT",
#     "AF007_SB_SC1_B": "UT",
#     "AF007_SB_SC1_C": "UT",
#     "AF007_SB_SC1_D": "INF",
#     "AF007_SB_SC1_E": "INF",
#     "AF007_SB_SC1_F": "INF",
# }
# Each sub-sample as one pseudo-bulk vector (3 UT + 3 INF replicates)
# sb_pb_df, sb_meta_df = process_pseudobulk(sb_samples)
# sb_pb_df.to_csv("pseudobulk_SB_gene_ids.csv")
# sb_meta_df.to_csv("pseudobulk_SB_metadata_new.csv")


def process_pseudobulk(samples):

    """
    Generate pseudobulk counts and metadata.
    Collect gene symbol and Ensembl ID mappings from all samples

    Parameters
    ----------
    samples : dict
        Sample names mapped to experimental conditions, e.g.
        {
            "AF007_SB_SC1_A": "UT",
            "AF007_SB_SC1_B": "UT",
            "AF007_SB_SC1_B": "UT",
            "AF007_SB_SC1_D": "INF",
            "AF007_SB_SC1_E": "INF",
            "AF007_SB_SC1_F": "INF"
        }

    Returns
    -------
    tuple
        Pseudobulk counts DataFrame and sample metadata DataFrame.
    """
    columns = {}
    meta_rows = []
    gene_maps = []

    for sample_name, condition in samples.items():
        path = h5ad_dir / f"{sample_name}_processed.h5ad"
        print(f"Reading {path.name} ...")
        ad = sc.read_h5ad(path)
        counts_sum = np.asarray(ad.layers["counts"].sum(axis=0)).flatten()
        columns[sample_name] = pd.Series(counts_sum, index=ad.var_names, name=sample_name)
        try:
            # Collect Ensembl IDs
            gene_maps.append(pd.Series(ad.var["gene_ids"].values, index=ad.var_names))
        except KeyError:
            gene_maps.append(pd.Series(ad.var["gene_id"].values, index=ad.var_names))

        meta_rows.append({"sample": sample_name, "condition": condition})
        print(f"  {sample_name}: {ad.n_obs} cells, {ad.n_vars} genes")

        # Combine counts
    sb_pb_df = pd.DataFrame(columns).fillna(0).astype(int)

    # Add gene_id column
    gene_ids = pd.concat(gene_maps)
    gene_ids = gene_ids[~gene_ids.index.duplicated(keep="first")]

    sb_pb_df.insert(0, "gene_id", sb_pb_df.index.map(gene_ids))
    sb_meta_df = pd.DataFrame(meta_rows).set_index("sample")

    print("returning pseudobulk data and metadata")

    return sb_pb_df, sb_meta_df



