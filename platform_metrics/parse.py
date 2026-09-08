import scanpy as sc
import pandas as pd
from platform_metrics.common import calculate_qc_metrics, get_stats_dict, get_cell_stats_df
from loguru import logger


def load_pa_sample(data_dir, project, sample):
    """
    Parse specific adata loading
    """

    adata = sc.read_mtx(
        data_dir / project / sample / "count" / "count_matrix.mtx.gz"
    )

    # reading in gene and cell data
    gene_data = pd.read_csv(
        data_dir / project / sample / "count" / "all_genes.csv.gz"
    )
    cell_meta = pd.read_csv(
        data_dir / project / sample / "count" / "cell_metadata.csv.gz"
    )

    ## Add gene metadata
    gene_data_use = gene_data.copy()

    gene_data_use["gene_name"] = gene_data_use["gene_name"].astype(str)
    gene_data_use = gene_data_use.set_index("gene_name")
    gene_data_use.index = gene_data_use.index.astype(str)
    gene_data_use.index.name = None

    adata.var = gene_data_use

    ## Add cell metadata
    cell_meta_cp = cell_meta.copy()

    cell_meta_cp["bc_wells"] = cell_meta_cp["bc_wells"].astype(str)
    cell_meta_cp = cell_meta_cp.set_index("bc_wells")
    cell_meta_cp.index = cell_meta_cp.index.astype(str)
    cell_meta_cp.index.name = None

    adata.obs = cell_meta_cp
    # Add the sample name to the adata object

    adata.obs["sample_name"] = sample

    logger.info(f'Returning adata for {sample}')

    return adata


def platform_metrics_df_pa(project_sample_dic, data_dir):
    """
    Parameters
    ----------
    project_sample_dic: Project and sample structure for file
    data_dir: Data directory

    Returns
    -------
    DF of parameters found in the adata files
    """

    metric_list = []
    for project, samples in project_sample_dic.items():
        for sample in samples:

            adata = load_pa_sample(data_dir, project, sample)
            adata = calculate_qc_metrics(adata)
            metric_dic = get_stats_dict(adata, sample, "Parse")
            metric_list.append(metric_dic)

    metric_df = pd.DataFrame(metric_list)
    logger.info(f'Returning {metric_df}')

    return metric_df

def get_cell_metrics_df_pa(project_sample_dic, data_dir, processed=False):
    """
    Parameters
    ----------
    project_sample_dic: Project and sample structure for file
    data_dir: Data directory

    Returns
    -------
    DF of parameters found in the adata files
    """

    all_cell_stats = []

    for project, samples in project_sample_dic.items():
        for sample in samples:

            if processed:
                adata = sc.read_h5ad(
                    data_dir / f"{sample}_raw.h5ad"
                )
            else:
                adata = load_pa_sample(data_dir, project, sample)
                adata = calculate_qc_metrics(adata)

            cell_stats = get_cell_stats_df(
                adata,
                sample_name=sample,
                platform="Parse"
            )

            all_cell_stats.append(cell_stats)

    cell_stats_df = pd.concat(all_cell_stats, ignore_index=True)

    logger.info(f"Returning cell stats dataframe with {len(cell_stats_df)} cells")

    return cell_stats_df

