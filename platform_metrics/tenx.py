import scanpy as sc
import pandas as pd
from platform_metrics.common import  get_stats_dict, get_cell_stats_df, calculate_qc_metrics
from loguru import logger


def load_10x_sample(data_dir, project, sample):
    adata = sc.read_10x_mtx(
        data_dir / project / sample / "count" / "sample_filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=True
    )
    # Add the sample name to the adata object
    adata.obs["sample_name"] = sample

    return adata


# def calculate_10x_qc_metrics(adata):
#     """
#     Small method to add qc metrics to tenx data
#     Param: adata
#     Returns: The adata object with the added qc metrics
#     """
#     adata.var_names_make_unique(join="_dup_")
#     sc.pp.calculate_qc_metrics(adata, inplace=True)
#
#     return adata

def platform_metrics_df_10x(project_sample_dic, data_dir):
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

            adata = load_10x_sample(data_dir, project, sample)
            adata = calculate_qc_metrics(adata)
            metric_dic = get_stats_dict(adata, sample, "10x")
            metric_list.append(metric_dic)

    metric_df = pd.DataFrame(metric_list)
    logger.info(f'Returning {metric_df}')

    return metric_df


def get_cell_metrics_df_10x(project_sample_dic, data_dir, processed=False):
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
                adata = load_10x_sample(data_dir, project, sample)
                adata = calculate_qc_metrics(adata)

            cell_stats = get_cell_stats_df(
                adata,
                sample_name=sample,
                platform="10x"
            )

            all_cell_stats.append(cell_stats)

    cell_stats_df = pd.concat(all_cell_stats, ignore_index=True)

    logger.info(f"Returning cell stats dataframe with {len(cell_stats_df)} cells")

    return cell_stats_df