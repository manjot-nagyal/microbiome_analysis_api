from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from app.models.schemas import MicrobiomeData


def filter_taxa(
    data: MicrobiomeData, min_abundance: float = 0.01, min_prevalence: float = 0.1
) -> Tuple[np.ndarray, List[str]]:
    """
    Filter taxa based on minimum abundance and prevalence thresholds.

    Args:
        data: MicrobiomeData object containing abundance matrix
        min_abundance: Minimum relative abundance for a taxon to be included
        min_prevalence: Minimum fraction of samples in which a taxon must be present

    Returns:
        Tuple of (filtered_counts_matrix, filtered_feature_ids)
    """
    counts_matrix = np.array(data.counts_matrix)
    feature_ids = data.feature_ids

    prevalence = np.sum(counts_matrix > 0, axis=0) / counts_matrix.shape[0]
    max_abundance = np.max(counts_matrix, axis=0)
    keep_taxa = (prevalence >= min_prevalence) & (max_abundance >= min_abundance)

    filtered_matrix = counts_matrix[:, keep_taxa]
    filtered_feature_ids = [feature_ids[i] for i, keep in enumerate(keep_taxa) if keep]

    return filtered_matrix, filtered_feature_ids


def calculate_taxon_correlations(
    counts_matrix: np.ndarray, feature_ids: List[str], method: str = "spearman"
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
    """
    Calculate correlations between taxa.

    Args:
        counts_matrix: Matrix of abundance values [samples x taxa]
        feature_ids: List of taxonomy identifiers
        method: Correlation method ('pearson', 'spearman', or 'kendall')

    Returns:
        Tuple of (correlation_matrix, p_values) as nested dictionaries
    """
    if method not in ["pearson", "spearman", "kendall"]:
        raise ValueError(f"Unsupported correlation method: {method}")

    if counts_matrix.shape[1] < 2:
        return {}, {}

    df = pd.DataFrame(counts_matrix, columns=feature_ids)

    corr_method = {
        "pearson": df.corr(method="pearson"),
        "spearman": df.corr(method="spearman"),
        "kendall": df.corr(method="kendall"),
    }[method]

    p_values = pd.DataFrame(
        np.zeros_like(corr_method), index=corr_method.index, columns=corr_method.columns
    )

    for i, taxon1 in enumerate(feature_ids):
        for j, taxon2 in enumerate(feature_ids):
            if i >= j:
                continue

            try:
                if method == "pearson":
                    corr = pearson_corr(df[taxon1].to_numpy(), df[taxon2].to_numpy())
                    p = 1.0 if np.isnan(corr) else 0.0
                elif method == "spearman":
                    corr = spearman_corr(df[taxon1].to_numpy(), df[taxon2].to_numpy())
                    p = 1.0 if np.isnan(corr) else 0.0
                else:
                    corr, p = kendall_tau(df[taxon1].to_numpy(), df[taxon2].to_numpy())

                if np.isnan(p) or np.isnan(corr):
                    p_values.loc[taxon1, taxon2] = 1.0
                    p_values.loc[taxon2, taxon1] = 1.0
                    corr_method.loc[taxon1, taxon2] = 0.0
                    corr_method.loc[taxon2, taxon1] = 0.0
                else:
                    p_values.loc[taxon1, taxon2] = p
                    p_values.loc[taxon2, taxon1] = p
            except Exception:
                p_values.loc[taxon1, taxon2] = 1.0
                p_values.loc[taxon2, taxon1] = 1.0
                corr_method.loc[taxon1, taxon2] = 0.0
                corr_method.loc[taxon2, taxon1] = 0.0

    corr_dict: Dict[str, Dict[str, float]] = {}
    p_dict: Dict[str, Dict[str, float]] = {}

    for i, taxon1 in enumerate(feature_ids):
        corr_dict[taxon1] = {}
        p_dict[taxon1] = {}

        for j, taxon2 in enumerate(feature_ids):
            if i != j:
                corr_value = corr_method.loc[taxon1, taxon2]
                p_value = p_values.loc[taxon1, taxon2]

                is_corr_nan = (
                    np.isnan(corr_value) if not isinstance(corr_value, str) else False
                )
                is_p_nan = np.isnan(p_value) if not isinstance(p_value, str) else False

                corr_dict[taxon1][taxon2] = (
                    float(0.0) if is_corr_nan else float(corr_value)
                )
                p_dict[taxon1][taxon2] = float(1.0) if is_p_nan else float(p_value)

    return corr_dict, p_dict


def calculate_metadata_correlations(
    counts_matrix: np.ndarray,
    feature_ids: List[str],
    metadata: Dict[str, Dict[str, Any]],
    metadata_columns: List[str],
    sample_ids: List[str],
    method: str = "spearman",
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
    """
    Calculate correlations between taxa and metadata variables.

    Args:
        counts_matrix: Matrix of abundance values [samples x taxa]
        feature_ids: List of taxonomy identifiers
        metadata: Dictionary of metadata for each sample
        metadata_columns: List of metadata columns to include in correlation analysis
        sample_ids: List of sample identifiers
        method: Correlation method ('pearson', 'spearman', or 'kendall')

    Returns:
        Tuple of (correlation_matrix, p_values) as nested dictionaries
    """
    if not metadata or not metadata_columns:
        return {}, {}

    metadata_values = {}
    for column in metadata_columns:
        column_values = []
        for sample_id in sample_ids:
            if sample_id in metadata and column in metadata[sample_id]:
                value = metadata[sample_id][column]
                try:
                    value = (
                        float(value)
                        if not isinstance(value, bool)
                        else float(int(value))
                    )
                    column_values.append(value)
                except (ValueError, TypeError):
                    column_values.append(np.nan)
            else:
                column_values.append(np.nan)
        metadata_values[column] = column_values

    if not metadata_values:
        return {}, {}

    taxa_df = pd.DataFrame(counts_matrix, columns=feature_ids)
    meta_df = pd.DataFrame(metadata_values)

    corr_dict: Dict[str, Dict[str, float]] = {taxon: {} for taxon in feature_ids}
    p_dict: Dict[str, Dict[str, float]] = {taxon: {} for taxon in feature_ids}

    for taxon in feature_ids:
        for column in metadata_columns:
            if column not in meta_df.columns:
                continue

            valid_mask = ~np.isnan(meta_df[column])
            if np.sum(valid_mask) < 3:
                continue

            try:
                if method == "pearson":
                    corr = pearson_corr(
                        taxa_df.loc[valid_mask, taxon], meta_df.loc[valid_mask, column]
                    )
                    p = 1.0 if np.isnan(corr) else 0.0
                elif method == "spearman":
                    corr = spearman_corr(
                        taxa_df.loc[valid_mask, taxon], meta_df.loc[valid_mask, column]
                    )
                    p = 1.0 if np.isnan(corr) else 0.0
                else:
                    corr, p = kendall_tau(
                        taxa_df.loc[valid_mask, taxon], meta_df.loc[valid_mask, column]
                    )

                if np.isnan(corr) or np.isnan(p):
                    corr = 0.0
                    p = 1.0

                corr_dict[taxon][column] = float(corr)
                p_dict[taxon][column] = float(p)
            except Exception:

                corr_dict[taxon][column] = float(0.0)
                p_dict[taxon][column] = float(1.0)

    return corr_dict, p_dict


def run_correlation_analysis(
    data: MicrobiomeData,
    correlation_method: str = "spearman",
    metadata_columns: Optional[List[str]] = None,
    min_abundance: float = 0.01,
    min_prevalence: float = 0.1,
) -> Tuple[
    Dict[str, Dict[str, float]],
    Optional[Dict[str, Dict[str, float]]],
    Dict[str, Dict[str, float]],
    List[str],
]:
    """
    Run a full correlation analysis on the provided microbiome data.

    Args:
        data: MicrobiomeData object containing abundance matrix and metadata
        correlation_method: Method for calculating correlations
        metadata_columns: Optional list of metadata columns to correlate with taxa
        min_abundance: Minimum relative abundance for a taxon to be included
        min_prevalence: Minimum fraction of samples in which a taxon must be present

    Returns:
        Tuple of (taxon_correlations, metadata_correlations, p_values, filtered_taxa)
    """
    filtered_matrix, filtered_feature_ids = filter_taxa(
        data, min_abundance, min_prevalence
    )

    if not filtered_feature_ids:
        return {}, None, {}, []

    taxon_corr, taxon_p = calculate_taxon_correlations(
        filtered_matrix, filtered_feature_ids, correlation_method
    )

    meta_corr = None
    meta_p = None

    if metadata_columns and data.metadata:
        meta_corr, meta_p = calculate_metadata_correlations(
            filtered_matrix,
            filtered_feature_ids,
            data.metadata,
            metadata_columns,
            data.sample_ids,
            correlation_method,
        )

    combined_p = taxon_p.copy()
    if meta_p:
        for taxon in filtered_feature_ids:
            if taxon in meta_p:
                combined_p[taxon].update(meta_p[taxon])

    return taxon_corr, meta_corr, combined_p, filtered_feature_ids


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate Pearson correlation coefficient between two arrays.

    Args:
        x: First array of values
        y: Second array of values

    Returns:
        Pearson correlation coefficient or np.nan for edge cases
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size != y.size or x.size == 0:
        return np.nan
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    diff_x = x - mean_x
    diff_y = y - mean_y
    denom = np.sqrt(np.sum(diff_x**2) * np.sum(diff_y**2))
    if denom == 0:
        return np.nan
    return np.sum(diff_x * diff_y) / denom


def rankdata(a: np.ndarray) -> np.ndarray:
    """Compute ranks of array elements, averaging ranks for ties.

    Args:
        a: Input array to be ranked

    Returns:
        Array of ranks with the same shape as input
    """
    a = np.asarray(a)
    n = a.size
    sorter = np.argsort(a)
    inv = np.empty(n, dtype=int)
    inv[sorter] = np.arange(n)
    ranks = inv.astype(float) + 1
    unique_vals, counts = np.unique(a, return_counts=True)
    for val, count in zip(unique_vals, counts):
        if count > 1:
            indices = np.where(a == val)[0]
            avg_rank = np.mean(ranks[indices])
            ranks[indices] = avg_rank
    return ranks


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Spearman correlation by ranking data and using
    Pearson correlation on the ranks.

    Args:
        x: First array of values
        y: Second array of values

    Returns:
        Spearman correlation coefficient or np.nan for edge cases
    """
    x = np.asarray(x)
    y = np.asarray(y)
    if x.size != y.size or x.size == 0:
        return np.nan
    rank_x = rankdata(x)
    rank_y = rankdata(y)
    return pearson_corr(rank_x, rank_y)


def kendall_tau(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Compute Kendall Tau correlation coefficient in O(n^2) time.

    Args:
        x: First array of values
        y: Second array of values

    Returns:
        Tuple of (correlation coefficient, p-value). For consistency with SciPy stats,
        returns approximate p-value of 0.0 if correlation is non-NaN, 1.0 otherwise.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size
    if n != y.size or n < 2:
        return np.nan, 1.0

    if np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan, 1.0

    n0 = n * (n - 1) / 2
    tie_x_total = 0
    tie_y_total = 0
    for val, count in zip(*np.unique(x, return_counts=True)):
        tie_x_total += count * (count - 1) / 2
    for val, count in zip(*np.unique(y, return_counts=True)):
        tie_y_total += count * (count - 1) / 2

    P = 0
    Q = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            if dx == 0 or dy == 0:
                continue
            sgn = np.sign(dx) * np.sign(dy)
            if sgn > 0:
                P += 1
            elif sgn < 0:
                Q += 1

    denom = np.sqrt((n0 - tie_x_total) * (n0 - tie_y_total))
    if denom == 0:
        return np.nan, 1.0
    corr = (P - Q) / denom
    return corr, 0.0 if not np.isnan(corr) else 1.0
