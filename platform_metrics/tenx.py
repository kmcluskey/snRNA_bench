import scanpy as sc
import pandas as pd
from platform_metrics.common import  get_stats_dict
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


def calculate_10x_qc_metrics(adata):
    """
    Small method to add qc metrics to tenx data
    Param: adata
    Returns: The adata object with the added qc metrics
    """
    adata.var_names_make_unique(join="_dup_")
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    return adata

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
            adata = calculate_10x_qc_metrics(adata)
            metric_dic = get_stats_dict(adata, sample, "10x")
            metric_list.append(metric_dic)

    metric_df = pd.DataFrame(metric_list)
    logger.info(f'Returning {metric_df}')

    return metric_df
