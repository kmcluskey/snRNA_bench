import scanpy as sc
import pandas as pd
from loguru import logger
from platform_metrics.common import calculate_qc_metrics, get_stats_dict


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
            metric_dic = get_stats_dict(adata, sample, "Scale")
            metric_list.append(metric_dic)

    metric_df = pd.DataFrame(metric_list)
    logger.info(f'Returning {metric_df}')

    return metric_df
