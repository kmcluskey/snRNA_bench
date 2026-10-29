from pathlib import Path

SERVER_ROOT = Path("/Volumes/kml29j/snRNA_bench")
pseudo_bulk_dir = SERVER_ROOT / "data" / "PseudoBulk" / "run_with_gene_ids"

## This was run on Pseudo_bulk-dseq jupyter notebook after the pseudobulk was run. Quan gave the
## actual bulk analysis.

# !/usr/bin/env python3
# deseq2.py

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


platforms = ["10X", "PA1", "SB"]


def filter_genes(counts_df, platform):
    n_before = counts_df.shape[1]

    # 1. remove 7SK
    mask_7sk = ~counts_df.columns.str.startswith("7SK")
    counts_df = counts_df.loc[:, mask_7sk]

    # 2. Remove low count gene
    min_samples = counts_df.shape[0] / 2
    mask_expr = (counts_df > 10).sum(axis=0) >= min_samples
    counts_df = counts_df.loc[:, mask_expr]

    n_after = counts_df.shape[1]
    print(f"  [{platform}] Genes before: {n_before} → after: {n_after} (removed {n_before - n_after})")
    return counts_df


def run_deseq2():

    for platform in platforms:

        print(f"\n{'=' * 50}")
        print(f"Platform: {platform}")
        print(f"{'=' * 50}")

        pb = pd.read_csv(
            pseudo_bulk_dir / f"pseudobulk_{platform}_gene_ids.csv",
            index_col=0
        )

        gene_ids = pb.pop("gene_id")
        counts = pb.T
        metadata = pd.read_csv(
            pseudo_bulk_dir / f"pseudobulk_{platform}_metadata_new.csv",
            index_col=0
        )
        metadata = metadata.loc[counts.index]
        metadata["condition"] = pd.Categorical(metadata["condition"], categories=["UT", "INF"])

        counts = filter_genes(counts, platform)

        # Replace gene symbols with Ensembl IDs
        counts.columns = gene_ids.loc[counts.columns].values

        if counts.columns.isna().any():
            raise ValueError(f"{platform}: Missing Ensembl IDs")

        if counts.columns.duplicated().any():
            raise ValueError(f"{platform}: Duplicate Ensembl IDs")

        print(f"  Counts shape: {counts.shape}")
        print(metadata)

        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design="~ condition",
            refit_cooks=True,
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=["condition", "INF", "UT"])
        stat_res.summary()

        results = stat_res.results_df.sort_values("padj")
        # Restore original gene symbols alongside Ensembl IDs
        ens2sym = pd.Series(gene_ids.index, index=gene_ids.values)

        results.index.name = "gene_id"
        results.insert(0, "symbol", results.index.map(ens2sym))

        print(results)

        out = f"DEA_{platform}_INF_vs_UT.csv"
        results.to_csv(out)

        print(f"Saved → {out}")


if __name__ == "__main__":
    run_deseq2()