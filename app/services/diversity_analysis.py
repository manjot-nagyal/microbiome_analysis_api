from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from scipy.spatial import distance

from app.models.schemas import MicrobiomeData


def shannon_diversity(counts: np.ndarray) -> float:
    """

    Calculate Shannon diversity index (H').
    Shannon diversity index takes into account both abundance
    and evenness of the species present.

    Formula: H' = -sum(p_i * ln(p_i)) for i = 1 to n
    where p_i is the proportion of individuals in the ith species.

    Args:
        data (np.ndarray): A 2D array where each row represents an
        individual and each column represents a species.

    Returns:
        Shannon diversity index value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    """
    counts = np.asarray(counts, dtype=float)

    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")

    counts = counts[counts > 0]

    if counts.size == 0:
        return np.nan

    proportions = counts / counts.sum()
    return -np.sum(proportions * np.log(proportions))


def simpson_diversity(counts: np.ndarray) -> float:
    """
    Calculate Simpson diversity index (1-D).

    The Simpson diversity index represents the probability that two randomly
    selected individuals in the habitat belong to different species. It ranges
    from 0 (no diversity) to 1 (infinite diversity).

    Formula: 1 - sum(p_i^2) where p_i is the proportion of species i.

    Args:
        counts: Array of abundance values for each taxon in a sample

    Returns:
        Simpson diversity index value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    """
    counts = np.asarray(counts, dtype=float)

    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")

    counts = counts[counts > 0]
    proportions = counts / counts.sum()

    if counts.size == 0:
        return np.nan

    proportions = counts / counts.sum()
    return float(1.0 - np.sum(proportions**2))


def pielou_evenness(counts: np.ndarray) -> float:
    """
    Pielou's evenness measures how evenly the species are distributed
    in a sample.
    It ranges from 0 (no evenness i.e. no dominance by one species) to
    1 (perfect evenness).

    Formula: J' = H' / ln(S) where H' is the Shannon diversity index and
    S is the number of species.

    Args:
        counts: Array of abundance values for each taxon in a sample

    Returns:
        Pielou's evenness value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    """
    counts = np.asarray(counts, dtype=float)

    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")

    if counts.size == 0:
        return np.nan

    total = counts.sum()
    if total == 0:
        return np.nan

    observed_species = np.count_nonzero(counts)

    if observed_species == 1:
        return 1.0

    proportions = counts / total
    # ignore proportions that are too tiny to be considered
    tiny = np.finfo(float).tiny
    pos_proortions = proportions[proportions > tiny]
    # calculate shannon diversity index
    h_prime = -np.sum(pos_proortions * np.log(pos_proortions))

    if np.isnan(h_prime) or np.isinf(h_prime):
        return np.nan

    max_diversity = np.log(observed_species) if observed_species > 1 else 0.0

    if h_prime == 0.0 or max_diversity == 0.0:
        return 0.0

    return h_prime / max_diversity


def chao1_estimator(counts: np.ndarray) -> float:
    """
    Calculate Chao1 estimator for species richness.

    Chao1 score is a non-parametric estimation of the number of species
    in a sample.
    It is based on the number of observed species and the number of rare
    species.

    Formula: Chao1 = S + [(f1(f1 - 1)) / (2(f2 + 1))]
    where S is the number of observed species, f1 are singletons and f2
    are doubletons.

    Args:
        counts: Array of abundance values for each taxon in a sample

    Returns:
        Chao1 estimator value, or NaN for edge cases

    Raises:
        ValueError: If any count value is negative
    """
    counts = np.asarray(counts, dtype=float)

    if np.any(counts < 0):
        raise ValueError("Negative abundance values are not allowed")
    if counts.size == 0:
        return np.nan
    counts = counts[counts > 0]
    S = np.count_nonzero(counts)
    if S == 0:
        return 0.0

    f1 = np.count_nonzero(counts == 1)
    f2 = np.count_nonzero(counts == 2)

    if f1 == 0:
        return float(S)
    if f2 == 0:
        return float(S) + (f1 * (f1 - 1)) / 2.0

    return float(S) + (f1 * (f1 - 1)) / (2.0 * (f2 + 1))


def calculate_alpha_diversity(
    data: MicrobiomeData,
    metrics: List[str] = ["shannon", "simpson", "pielou", "chao1"],
) -> Dict[str, float]:
    """
    Calculate alpha diversity metrics for each sample

    Args:
        data: MicrobiomeData object containing abundance matrix
        metrics: List of alpha diversity metrics to calculate

    Returns:
        Dictionary of alpha diversity metrics for each sample

    Raises:
        ValueError: If any unsupported metric is provided
    """

    counts_matrix = np.array(data.counts_matrix)
    sample_ids = data.sample_ids

    results = {}

    for metric in metrics:
        metric_results = {}

        for j, sample_id in enumerate(sample_ids):
            sample_counts = counts_matrix[j, :]
            if metric == "shannon":
                alpha_diversity_value = shannon_diversity(sample_counts)
            elif metric == "simpson":
                alpha_diversity_value = simpson_diversity(sample_counts)
            elif metric == "pielou":
                alpha_diversity_value = pielou_evenness(sample_counts)
            elif metric == "chao1":
                alpha_diversity_value = chao1_estimator(sample_counts)
            else:
                raise ValueError(f"Unsupported metric: {metric}")

            metric_results[sample_id] = alpha_diversity_value

        results[metric] = metric_results

    return results


def calculate_beta_diversity(
    data: MicrobiomeData,
    metrics: List[str] = [
        "braycurtis",
        "jaccard",
    ],
) -> Dict[str, List[List[float]]]:
    """
    Calculate beta diversity distance matrices between samples

    Args:
        data: MicrobiomeData object containing abundance matrix
        metrics: List of beta diversity metrics to calculate

    Returns:
        Dictionary of beta diversity distance matrices for each sample

    Raises:
        ValueError: If any unsupported metric is provided
    """

    counts_matrix = np.array(data.counts_matrix)
    results = {}

    for metric in metrics:
        if metric == "braycurtis":
            dist = distance.pdist(counts_matrix, metric="braycurtis")
            results[metric] = distance.squareform(dist).tolist()
        elif metric == "jaccard":
            binary_matrix = (counts_matrix > 0).astype(float)
            dist = distance.pdist(binary_matrix, metric="jaccard")
            results[metric] = distance.squareform(dist).tolist()
        else:
            raise ValueError(f"Unsupported metric: {metric}")

    return results


def calculate_group_comparison(
    data: MicrobiomeData,
    group_by: str,
    diversity_metrics: Dict[str, Dict[str, float]],
) -> Optional[Dict[str, Dict[str, Dict[str, float]]]]:
    """
    Compare alpha diversity metrics between groups using statistical tests.

    This function performs pairwise comparisons between sample groups
    using Mann-Whitney U test. This will determine if diversity differs
    significantly between groups.

    Args:
        data: MicrobiomeData object containing abundance matrix and metadata
        group_by: Metadata column to group samples by for comparison
        diversity_metrics: Dictionary of alpha diversity metrics for each sample
    Returns:
        Dictionary of group comparison results for each alpha diversity metric
    """

    if not data.metadata or not group_by:
        return None

    groups: Dict[str, List[str]] = {}

    for sample_id in data.sample_ids:
        if sample_id in data.metadata and group_by in data.metadata[sample_id]:
            group_value = data.metadata[sample_id][group_by]
            if group_value not in groups:
                groups[group_value] = []
            groups[group_value].append(sample_id)

    if len(groups) < 2:
        return None

    group_comparison = {}

    for metric, values in diversity_metrics.items():
        metric_results = {}

        group_keys = list(groups.keys())
        for i in range(len(group_keys)):
            for j in range(i + 1, len(group_keys)):
                group_a = str(group_keys[i])
                group_b = str(group_keys[j])
                group_a_values = [
                    values[sample_id]
                    for sample_id in groups[group_a]
                    if sample_id in values
                ]
                group_b_values = [
                    values[sample_id]
                    for sample_id in groups[group_b]
                    if sample_id in values
                ]

                if len(group_a_values) < 1 or len(group_b_values) < 1:
                    continue

                stat, p_value = stats.mannwhitneyu(group_a_values, group_b_values)

                comparison_key = f"{group_a}_vs_{group_b}"

                metric_results[comparison_key] = {
                    "p_value": float(p_value),
                    "statistic": float(stat),
                    "significant": p_value < 0.05,
                    f"{group_a}_mean": np.mean(group_a_values),
                    f"{group_b}_mean": np.mean(group_b_values),
                }

        group_comparison[metric] = metric_results

    return group_comparison


def run_diversity_analysis(
    data: MicrobiomeData,
    metrics: List[str] = ["shannon", "simpson", "pielou", "chao1"],
    group_by: Optional[str] = None,
) -> Tuple[
    Dict[str, Dict[str, float]],
    Dict[str, List[List[float]]],
    Optional[Dict[str, Dict[str, Dict[str, float]]]],
]:
    """
    Run a full diversity analysis on the provided microbiome data.

    This function orchestrates the complete diversity analysis workflow by:
    1. Calculating alpha diversity metrics for each sample
    2. Computing beta diversity distance matrices between samples
    3. Performing statistical comparisons between groups if metadata
    is provided

    Args:
        data: MicrobiomeData object containing abundance matrix and metadata
        metrics: List of alpha diversity metrics to calculate
        group_by: Optional metadata column to group samples by for comparison

    Returns:
        Tuple of (alpha_diversity, beta_diversity, group_comparison) results
    """

    alpha_diversity = calculate_alpha_diversity(data, metrics)

    beta_diversity = calculate_beta_diversity(data)

    group_comparison = calculate_group_comparison(data, group_by, alpha_diversity)

    return alpha_diversity, beta_diversity, group_comparison
