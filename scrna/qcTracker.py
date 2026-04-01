import pandas as pd
from loguru import logger
import numpy as np

class QCTracker:


    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        self.records = []
        self.prev_cells = None
        self.prev_genes = None
        self.step_names = set()

    def _median_or_none(self, series):
        if series is None:
            return None
        if len(series) == 0:
            return None
        return float(np.nanmedian(series))


    def record(self, adata, step_name, param_dict=None):

        if step_name in self.step_names:
            print(f"Warning: step '{step_name}' already exists. No changes made.")
            return

        self.step_names.add(step_name)

        n_cells = adata.n_obs
        n_genes = adata.n_vars

        if self.prev_cells is None:
            cells_removed = 0
            genes_removed = 0
        else:
            cells_removed = self.prev_cells - n_cells
            genes_removed = self.prev_genes - n_genes

        # Keep your nice param formatting
        if param_dict:
            params = ", ".join(param_dict.keys())
            values = ", ".join(str(x) for x in param_dict.values())
            print(params, values)
        else:
            params = None
            values = None

        total_umis = (
            int(adata.obs["total_counts"].sum())
            if "total_counts" in adata.obs
            else None
        )

        median_umis = (
            adata.obs["total_counts"].median().round(2)
            if "total_counts" in adata.obs
            else None
        )

        median_genes_per_cell = (
            adata.obs["n_genes_by_counts"].median().round(2)
            if "n_genes_by_counts" in adata.obs
            else None
        )

        # Safely get QC metrics if present
        median_mt = (
            adata.obs["pct_counts_mt"].median().round(2)
            if "pct_counts_mt" in adata.obs
            else None
        )

        median_rp = (
            adata.obs["pct_counts_rp"].median().round(2)
            if "pct_counts_rp" in adata.obs
            else None
        )

        if "predicted_doublet" in adata.obs:
            print("here")
            rate = adata.obs["predicted_doublet"].mean() * 100
            print ("yo", rate)
            if not np.isnan(rate):
                doublet_rate = round(rate, 2)
            else:
                doublet_rate = None
        else:
            doublet_rate = None

        # --- Cell type / cluster summary (only if present) ---

        if "leiden" in adata.obs:
            n_clusters = adata.obs["leiden"].nunique()
        else:
            n_clusters = None

        if "cell_type" in adata.obs:
            celltype_counts_series = adata.obs["cell_type"].value_counts()

            # Names of annotated types
            cell_types = ", ".join(celltype_counts_series.index.tolist())

            # Compact string with counts in brackets
            celltype_counts = ", ".join(
                f"{ct} ({n})"
                for ct, n in celltype_counts_series.items()
            )
        else:
            cell_types = None
            celltype_counts = None

        self.records.append({
            "dataset": self.dataset_name,
            "step": step_name,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "total_umis": total_umis,
            "median_umis": median_umis,
            "median_genes_cell": median_genes_per_cell,
            "params": params,
            "values": values,
            "cells_removed": cells_removed,
            "genes_removed": genes_removed,
            "median_pct_mt": median_mt,
            "median_pct_rp": median_rp,
            "doublet_rate": doublet_rate,
            "clusters": n_clusters,
            "cell_types": cell_types,
            "celltype_counts": celltype_counts,
        })

        self.prev_cells = n_cells
        self.prev_genes = n_genes

    def to_dataframe(self):
        return pd.DataFrame(self.records)