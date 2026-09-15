import scanpy as sc
import pandas as pd
from loguru import logger
from platform_metrics.common import calculate_qc_metrics, get_stats_dict, get_cell_stats_df


def load_sb_sample(data_dir, project, sample):
    """
    Given the data dir, project and sample return the adata file for scale data
    """
    adata = sc.read_10x_mtx(
        data_dir / project / sample / "count",
        var_names="gene_symbols",
        cache=True
    )
    ### Extra step scale bio data
    adata.layers["counts"] = adata.X.astype("int64").copy()

    # Add the sample name to the adata object
    adata.obs["sample_name"] = sample

    return adata

def platform_metrics_df_sb(project_sample_dic, data_dir):
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

            adata = load_sb_sample(data_dir, project, sample)
            adata = calculate_qc_metrics(adata)
            metric_dic = get_stats_dict(adata, sample, "Scalebio")
            metric_list.append(metric_dic)

    metric_df = pd.DataFrame(metric_list)
    logger.info(f'Returning {metric_df}')

    return metric_df


def get_cell_metrics_df_sb(project_sample_dic, data_dir, processed=False):
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
                    data_dir / f"{sample}_processed.h5ad"
                )
            else:
                adata = load_sb_sample(data_dir, project, sample)
                adata = calculate_qc_metrics(adata)

            cell_stats = get_cell_stats_df(
                adata,
                sample_name=sample,
                platform="Scalebio"
            )

            all_cell_stats.append(cell_stats)

    cell_stats_df = pd.concat(all_cell_stats, ignore_index=True)

    logger.info(f"Returning cell stats dataframe with {len(cell_stats_df)} cells")

    return cell_stats_df