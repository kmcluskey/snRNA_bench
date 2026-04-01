
import numpy as np
import pandas as pd
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score



def print_cells_and_genes(adata):
    n_cells, n_genes = adata.shape
    print(f"Dataset contains {n_cells:,} cells and {n_genes:,} genes.")

def test_clustering_plots_plus_rank(
    adata,
    pcs_list,
    neighbors_list,
    res_list,
    flavor=None,
    n_top_genes=2000,
    random_state=0,
    n_rank_genes=5,
    make_plots=True,
    use_rep=None,
):
    """
    Parameter sweep for clustering optimisation.

    Uses HVGs for PCA (without subsetting genes),
    fixes randomness for reproducibility,
    and returns a summary table to help choose parameters.

    Parameters
    ----------
    adata : AnnData
    pcs_list : list[int]
    neighbors_list : list[int]
    res_list : list[float]
    flavor : str or None
        HVG method (e.g. "cell_ranger", "seurat_v3")
    n_top_genes : int
    random_state : int
    min_cells_per_cluster : int
        Flags parameter sets that produce tiny clusters
    n_rank_genes : int
        Number of marker genes to visualise
    make_plots : bool

    Returns
    -------
    pandas.DataFrame summary of parameter performance
    """

    # Work on copy so original adata is untouched

    # # HVG selection (stored in base.var["highly_variable"])
    # if flavor:
    #     sc.pp.highly_variable_genes(
    #         base,
    #         flavor=flavor,
    #         n_top_genes=n_top_genes
    #     )
    # else:
    #     sc.pp.highly_variable_genes(
    #         base,
    #         n_top_genes=n_top_genes
    #     )
    #
    # # PCA using HVGs only (without subsetting genes)
    # sc.pp.pca(
    #     base,
    #     use_highly_variable=True,
    #     random_state=random_state
    # )
    print ('deleted test')
    results = []

    for pcs in pcs_list:
        for neighbors in neighbors_list:
            for res in res_list:

                adata_base = adata.copy()

                # adata_base = base.copy()

                print(f"Testing: PCs={pcs}, neighbors={neighbors}, resolution={res}")

                ## Calculating neighbours with or without harmony
                if use_rep:
                    sc.pp.neighbors(
                        adata_base,
                        n_neighbors=neighbors,
                        use_rep=use_rep,
                        random_state=random_state
                    )
                else:
                    sc.pp.neighbors(
                        adata_base,
                        n_neighbors=neighbors,
                        n_pcs=pcs,
                        random_state=random_state
                    )

                # Leiden clustering
                sc.tl.leiden(
                    adata_base,
                    resolution=res,
                    key_added="leiden",
                    random_state=random_state
                )

                # UMAP
                sc.tl.umap(adata_base, random_state=random_state)

                # Cluster stats
                n_clusters = adata_base.obs["leiden"].nunique()

                # Silhouette score on PCA embedding
                sil = np.nan
                if n_clusters > 1:
                    rep = use_rep if use_rep else "X_pca"
                    X = adata_base.obsm[rep][:, :pcs]
                    labels = adata_base.obs["leiden"].astype(str).to_numpy()
                    sil = float(silhouette_score(X, labels))

                # Rank genes
                sc.tl.rank_genes_groups(
                    adata_base,
                    "leiden",
                    method="wilcoxon"
                )

                if make_plots:
                    sc.pl.umap(
                        adata_base,
                        color="leiden",
                        title=f"PCs={pcs}, nn={neighbors}, res={res} | "
                              f"clusters={n_clusters} | sil={sil:.3f}"
                    )

                    sc.pl.rank_genes_groups_dotplot(
                    adata_base, n_genes=n_rank_genes, standard_scale="var",)   # optional but usually nicer

                results.append({
                    "pcs": pcs,
                    "neighbors": neighbors,
                    "resolution": res,
                    "n_clusters": n_clusters,
                    "silhouette_pca": sil
                })

    summary = pd.DataFrame(results)

    summary_table = summary.sort_values(
        by=["silhouette_pca", "n_clusters"],
        ascending=[False, True])

    return summary_table


def get_stats_table(adata, processed=False, title=None):

    """
    Just a wee table to compare the stats and produce a table for adata
    """
    summary_counts = adata.obs["total_counts"].describe()
    summary_genes = adata.obs["n_genes_by_counts"].describe()

    metrics = [
        "Number of cells",
        "Number of genes",
        "Median UMIs per cell",
        "Max UMIs per cell",
        "Median genes per cell",

    ]

    values = [
        adata.n_obs,
        adata.n_vars,
        int(summary_counts["50%"]),
        int(summary_counts["max"]),
        int(summary_genes["50%"])
    ]

    if processed:
        metrics.extend([
            "% Ribosomal Protein Counts",
            "% Mitochondrial Protein Counts",
            "% Predicted doublets"
        ])

        values.extend([
            adata.obs["pct_counts_rp"].mean().round(2),
            adata.obs["pct_counts_mt"].mean().round(2),
            adata.obs["predicted_doublet"].mean().round(2) * 100
        ])

    table = pd.DataFrame({
        "Metric": metrics,
        "Value": values
    })

    if title:
        table.attrs["title"] = title


    return table


def get_sig_ranked_df(df, lfc=1, score=10, pval=0.01):
    """
    params: df with headings logfoldchanges, scores and pvals_adj
    from sc.get.rank_genes_groups_df
    Return: df - but filtered to the conditions passed in

    """
    print(f'getting ranked df with lfc = {lfc}, score ={score} and pval = {pval}')
    sig_df = df[
        (df["logfoldchanges"] > lfc) &
        (df["scores"] > score) &
        (df["pvals_adj"] < pval)
        ]

    return sig_df

def get_gene_matches(name, gene_df, genes_to_match):
    """
    Take a DF with names and match to a list of genes - calculate the mean values for the several features
    and return the matching DF and a summary.
    Used for matching to list of interferon genes for example.
    """
    mean_lfc_comp  = gene_df["logfoldchanges"].mean()
    mean_scores_comp  = gene_df["scores"].mean()

    overlap = gene_df[gene_df.names.isin(genes_to_match)]

    n = len(overlap)
    mean_lfc = overlap["logfoldchanges"].mean()
    mean_score =  overlap["scores"].mean()

    wee_summary = pd.DataFrame([{
        "n_sig_genes_comp": gene_df.shape[0], #Sig genes in the comparison passes/entire DF
        "mean_lfc_comp": mean_lfc_comp,
        "mean_score_comp": mean_scores_comp,
        "comparison": name,
        "n_match_genes": n, # Matching genes between those passed and those in the dataset
        "mean_lfc_match": mean_lfc,
        "mean_score_match": mean_score
    }])

    return overlap, wee_summary