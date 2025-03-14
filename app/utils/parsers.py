import json
from io import StringIO
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from app.models.schemas import MicrobiomeData


def parse_csv_data(content_str: str) -> MicrobiomeData:
    """
    Parse the CSV data into a MicrobiomeData object
    """
    try:
        try:
            df = pd.read_csv(
                StringIO(content_str), index_col=0, na_values=["", "NA", "nan"]
            )

            df = df.dropna(axis=1, how="all")

            # Separate metadata from data
            metadata_cols = [col for col in df.columns if col.startswith("metadata_")]
            data_cols = [col for col in df.columns if not col.startswith("metadata_")]
            sample_ids = df.index.tolist()

            # Extract feature counts
            counts_df = df[data_cols]

            metadata: Optional[Dict[str, Dict[str, Any]]] = None
            if metadata_cols:
                metadata_df = df[metadata_cols]
                metadata = {}

                for sample_id in df.index:
                    metadata[sample_id] = {}
                    for meta_col in metadata_cols:
                        clean_meta_name = meta_col.replace("metadata_", "")
                        value = metadata_df.loc[sample_id, meta_col]
                        # Skip NaN values in metadata completely
                        if pd.isna(value):
                            continue
                        else:
                            metadata[sample_id][clean_meta_name] = value
            try:
                counts_matrix: List[List[float]] = []
                for sample in counts_df.values:
                    sample_values = []
                    for val in sample:

                        if pd.isna(val):
                            sample_values.append(0.0)
                        else:
                            sample_values.append(float(val))
                    counts_matrix.append(sample_values)
            except ValueError as e:
                raise ValueError("Error parsing CSV file: " + str(e))

            microbiome_data = MicrobiomeData(
                sample_ids=sample_ids,
                feature_ids=data_cols,
                counts_matrix=counts_matrix,
                metadata=metadata,
            )

            return microbiome_data
        except Exception as e:
            print(f"Error parsing CSV data: {str(e)}")
            raise
    except Exception as e:
        print(f"Error parsing CSV data: {str(e)}")
        raise


def parse_json_data(content: str) -> MicrobiomeData:
    """
    Parse JSON format microbiome data.

    Expected format:
    {
        "sample_ids": ["sample1", "sample2", ...],
        "feature_ids": ["taxon1", "taxon2", ...],
        "counts_matrix": [[val1, val2, ...], [val1, val2, ...], ...],
        "metadata": {
            "sample1": {"metadata1": value1, "metadata2": value2, ...},
            "sample2": {"metadata1": value1, "metadata2": value2, ...},
            ...
        }
    }

    Args:
        content: JSON data as string

    Returns:
        MicrobiomeData object with parsed data
    """
    data = json.loads(content)

    return MicrobiomeData(**data)


def normalize_counts(data: MicrobiomeData) -> MicrobiomeData:
    """
    Normalize counts data to relative abundance (sum to 1 for each sample).

    Args:
        data: MicrobiomeData object with raw abundance counts

    Returns:
        MicrobiomeData object with normalized abundance values
    """

    counts_matrix = np.array(data.counts_matrix)

    counts_matrix = np.nan_to_num(counts_matrix, nan=0.0)

    row_sums = counts_matrix.sum(axis=1, keepdims=True)

    row_sums[row_sums == 0] = 1.0
    normalized_matrix = counts_matrix / row_sums

    data_dict = data.model_dump()
    data_dict["counts_matrix"] = normalized_matrix.tolist()

    return MicrobiomeData(**data_dict)
