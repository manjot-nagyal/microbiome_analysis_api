# import numpy as np
import pandas as pd

from app.models.schemas import MicrobiomeData


def parse_csv_data(content_str: str) -> MicrobiomeData:
    """
    Parse the CSV data into a MicrobiomeData object
    """

    df = pd.read_csv(content_str, index_col=0)

    # Separate metadata from data
    metadata_cols = df.columns[df.columns.str.startswith("metadata_")]
    data_cols = df.columns[~df.columns.str.startswith("metadata_")]

    # Extract feature counts
    counts_df = df[data_cols]

    # TODO: Extract and format metadata
    md_df = df[metadata_cols]

    # Create MicrobiomeData object
    microbiome_data = MicrobiomeData(
        sample_ids=counts_df.index.tolist(),
        feature_ids=counts_df.columns.tolist(),
        counts_matrix=counts_df.values.tolist(),
        metadata=(md_df.to_dict(orient="records") if md_df is not None else None),
    )

    return microbiome_data
