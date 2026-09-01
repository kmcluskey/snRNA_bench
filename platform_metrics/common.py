
def get_stats_dict(adata, sample_name, platform):


    stats_dict = {
                "sample": sample_name,
                "platform": platform,
                "n_cells": adata.n_obs,
                "n_features": adata.n_vars,
                "median_genes_per_cell": adata.obs["n_genes_by_counts"].median(),
                "median_counts_per_cell": adata.obs["total_counts"].median(),
                "mean_genes_per_cell": adata.obs["n_genes_by_counts"].mean(),
                "mean_counts_per_cell": adata.obs["total_counts"].mean(),
            }

    print("returning", stats_dict)
