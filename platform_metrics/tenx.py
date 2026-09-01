import scanpy as sc


def load_10x_sample(data_dir, project, sample):
    adata = sc.read_10x_mtx(
        data_dir / project / sample / "count" / "sample_filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=True
    )

    return adata


def calculate_10x_qc_metrics(adata):
    """
    Small method to add qc metrics to tenx data
    Param: adata
    Returns: The adata object with the added qc metrics
    """
    adata.var_names_make_unique(join="_dup_")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    return adata
