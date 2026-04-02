import os
import re
import pandas as pd
from loguru import logger
from IPython.core.display_functions import display

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
UPSTREAM_DIR = PROJECT_ROOT / "data" / "upstream_data"

def read_upstream_csv(filename):

    csv_path = UPSTREAM_DIR / filename

    if not csv_path.exists():
        raise FileNotFoundError(f"File not found: {csv_path}")

    return pd.read_csv(csv_path)

def clean_upsteam_metrics(metrics_file_csv):
    """
    A class to take the processed output from 10x (cr), Parse (pb) and ScaleBio (sb)
    and put in a format useable for figures and table.
    """
    metrics_df = read_upstream_csv(metrics_file_csv)

    display(metrics_df)
