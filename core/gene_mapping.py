from loguru import logger
from pathlib import Path
import pandas as pd
import scanpy as sc

"""
A set of helper files that can be used to standardise the gene_mapping for human snRNA-seq and bulk_RNA-seq
"""


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = PROJECT_ROOT / "reference_data"

def get_ens2sym():

    """
    Functions for standardising gene identifiers across
    human snRNA-seq and bulk RNA-seq datasets.
    params: return: ens2sym a pandas series with the ens2sym relationship
    """


    raw = pd.read_csv(
        Path(REFERENCE_DIR) / "biomart_export.tsv",
        sep="\t",
        header=None,
        usecols=[2, 3],
        names=["ens", "sym"]
    )

    #keep the **first** occurrence of each Ensembl ID
    ens2sym = (
        raw[~raw["ens"].duplicated(keep="first")]   # drop later duplicates
            .set_index("ens")["sym"]
            .str.upper()
    )

    logger.info("Returning {} unique Ensembl IDs", len(ens2sym))

    return ens2sym


def get_gene_map(samples, h5ad_dir):

    '''
    ### Take all the single cell datasets from a platform in order to
    ### get a comprehensive gene map symbol:ENSG - can be used for one dataset.
    return: symbol: ENSG dataframe
    '''

    gene_maps = []

    for sample in samples:
        ad = sc.read_h5ad(h5ad_dir / f"AF007_{sample}_processed.h5ad", backed="r")

        col = "gene_ids" if "gene_ids" in ad.var.columns else "gene_id"

        gene_maps.append(pd.DataFrame({
            "symbol": ad.var_names.to_numpy(),
            "gene_id": ad.var[col].to_numpy()
        }))

        ad.file.close()

    gene_map = (
        pd.concat(gene_maps, ignore_index=True)
        .dropna(subset=["symbol", "gene_id"])
        .drop_duplicates()
    )

    print(f"Unique mappings: {len(gene_map):,}")

    # Check ambiguous symbols
    ambiguous = gene_map.groupby("symbol")["gene_id"].nunique()
    print(f"Ambiguous symbols: {(ambiguous > 1).sum()}")

    return gene_map


def clean_degs(deseg_df):
    """
    Standardise bulk and pseudobulk DESeq2 results.
    Retain original symbols and add BioMart symbols using Ensembl IDs.
    """
    deseg_df = deseg_df.copy()

    # Standardise column names (bulk RNA-seq)
    deseg_df = deseg_df.rename(columns={
        "Unnamed: 0": "gene_id",
        "geneName": "symbol"
    })

    # Remove Ensembl version numbers, if present
    deseg_df["gene_id"] = (
        deseg_df["gene_id"].astype("string")
        .str.replace(r"\.\d+$", "", regex=True)
    )

    # Standardise capitalisation
    deseg_df["symbol"] = deseg_df["symbol"].str.upper()

    # Add BioMart symbols
    ens2sym = get_ens2sym()
    deseg_df["biomart_symbol"] = deseg_df["gene_id"].map(ens2sym).str.upper()

    logger.info(
        "Cleaned DEGs: {} genes, {} missing symbols, "
        "{} missing Ensembl IDs, {} missing BioMart symbols",
        len(deseg_df),
        deseg_df["symbol"].isna().sum(),
        deseg_df["gene_id"].isna().sum(),
        deseg_df["biomart_symbol"].isna().sum(),
    )

    return deseg_df
