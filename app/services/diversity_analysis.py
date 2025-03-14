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
    try:
        counts = np.asarray(counts, dtype=float)

        print(f"Shannon diversity - Original counts: {counts[:5]}")
        print(f"Contains NaN: {np.isnan(counts).any()}")

        counts = np.nan_to_num(counts, nan=0.0)
        print(f"After nan_to_num: {counts[:5]}")

        if np.any(counts < 0):
            raise ValueError("Negative abundance values are not allowed")

        counts = counts[counts > 0]
        print(
            f"After filtering zeros: {counts[:5] if len(counts) > 0 else 'Empty array'}"
        )

        if counts.size == 0:
            print("Shannon diversity: Empty array after filtering, returning NaN")
            return np.nan

        proportions = counts / counts.sum()
        result = -np.sum(proportions * np.log(proportions))
        print(f"Shannon diversity result: {result}")
        return result
    except Exception as e:
        print(f"Error in shannon_diversity: {str(e)}")
        raise


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
    try:
        counts = np.asarray(counts, dtype=float)

        print(f"Simpson diversity - Original counts: {counts[:5]}")
        print(f"Contains NaN: {np.isnan(counts).any()}")

        counts = np.nan_to_num(counts, nan=0.0)
        print(f"After nan_to_num: {counts[:5]}")

        if np.any(counts < 0):
            raise ValueError("Negative abundance values are not allowed")

        counts = counts[counts > 0]
        print(
            f"After filtering zeros: {counts[:5] if len(counts) > 0 else 'Empty array'}"
        )

        if counts.size == 0:
            print("Simpson diversity: Empty array after filtering, returning NaN")
            return np.nan

        total = counts.sum()
        if total == 0:
            print("Simpson diversity: Total is zero, returning NaN")
            return np.nan

        proportions = counts / total
        result = float(1.0 - np.sum(proportions**2))
        print(f"Simpson diversity result: {result}")
        return result
    except Exception as e:
        print(f"Error in simpson_diversity: {str(e)}")
        raise


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
    try:
        counts = np.asarray(counts, dtype=float)

        print(f"Pielou evenness - Original counts: {counts[:5]}")
        print(f"Contains NaN: {np.isnan(counts).any()}")

        counts = np.nan_to_num(counts, nan=0.0)
        print(f"After nan_to_num: {counts[:5]}")

        if np.any(counts < 0):
            raise ValueError("Negative abundance values are not allowed")

        if counts.size == 0:
            print("Pielou evenness: Empty array, returning NaN")
            return np.nan

        total = counts.sum()
        if total == 0:
            print("Pielou evenness: Total is zero, returning NaN")
            return np.nan

        observed_species = np.count_nonzero(counts)
        print(f"Observed species: {observed_species}")

        if observed_species == 1:
            print("Pielou evenness: Only one species, returning 1.0")
            return 1.0

        proportions = counts / total
        # ignore proportions that are too tiny to be considered
        tiny = np.finfo(float).tiny
        pos_proortions = proportions[proportions > tiny]
        # calculate shannon diversity index
        h_prime = -np.sum(pos_proortions * np.log(pos_proortions))
        print(f"Shannon diversity (h_prime): {h_prime}")

        if np.isnan(h_prime) or np.isinf(h_prime):
            print("Pielou evenness: h_prime is NaN or Inf, returning NaN")
            return np.nan

        max_diversity = np.log(observed_species) if observed_species > 1 else 0.0
        print(f"Max diversity: {max_diversity}")

        if h_prime == 0.0 or max_diversity == 0.0:
            print("Pielou evenness: h_prime or max_diversity is zero, returning 0.0")
            return 0.0

        result = h_prime / max_diversity
        print(f"Pielou evenness result: {result}")
        return result
    except Exception as e:
        print(f"Error in pielou_evenness: {str(e)}")
        raise


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
    try:
        counts = np.asarray(counts, dtype=float)

        print(f"Chao1 estimator - Original counts: {counts[:5]}")
        print(f"Contains NaN: {np.isnan(counts).any()}")

        counts = np.nan_to_num(counts, nan=0.0)
        print(f"After nan_to_num: {counts[:5]}")

        if np.any(counts < 0):
            raise ValueError("Negative abundance values are not allowed")
        if counts.size == 0:
            print("Chao1 estimator: Empty array, returning NaN")
            return np.nan
        counts = counts[counts > 0]
        print(
            f"After filtering zeros: {counts[:5] if len(counts) > 0 else 'Empty array'}"
        )

        S = np.count_nonzero(counts)
        if S == 0:
            print("Chao1 estimator: S is zero, returning 0.0")
            return 0.0

        f1 = np.count_nonzero(counts == 1)
        f2 = np.count_nonzero(counts == 2)
        print(f"S: {S}, f1: {f1}, f2: {f2}")

        if f1 == 0:
            print("Chao1 estimator: f1 is zero, returning S: {S}")
            return float(S)
        if f2 == 0:
            result = float(S) + (f1 * (f1 - 1)) / 2.0
            print(f"Chao1 estimator (f2=0): {result}")
            return result

        result = float(S) + (f1 * (f1 - 1)) / (2.0 * (f2 + 1))
        print(f"Chao1 estimator result: {result}")
        return result
    except Exception as e:
        print(f"Error in chao1_estimator: {str(e)}")
        raise


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
    try:
        print("Starting alpha diversity calculation")
        counts_matrix = np.array(data.counts_matrix)
        sample_ids = data.sample_ids

        print(f"Counts matrix shape: {counts_matrix.shape}")
        print(f"Sample IDs: {sample_ids[:5]}...")

        results = {}

        for metric in metrics:
            print(f"Calculating {metric} diversity")
            metric_results = {}

            for j, sample_id in enumerate(sample_ids):
                print(f"Processing sample {j+1}/{len(sample_ids)}: {sample_id}")
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
            print(f"Completed {metric} diversity calculation")

        return results
    except Exception as e:
        print(f"Error in calculate_alpha_diversity: {str(e)}")
        raise


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
    try:
        print("Starting beta diversity calculation")
        counts_matrix = np.array(data.counts_matrix)

        print(f"Beta diversity - Contains NaN: {np.isnan(counts_matrix).any()}")

        counts_matrix = np.nan_to_num(counts_matrix, nan=0.0)

        results = {}

        for metric in metrics:
            print(f"Calculating {metric} beta diversity")
            if metric == "braycurtis":
                dist = distance.pdist(counts_matrix, metric="braycurtis")
                results[metric] = distance.squareform(dist).tolist()
            elif metric == "jaccard":
                binary_matrix = (counts_matrix > 0).astype(float)
                dist = distance.pdist(binary_matrix, metric="jaccard")
                results[metric] = distance.squareform(dist).tolist()
            else:
                raise ValueError(f"Unsupported metric: {metric}")
            print(f"Completed {metric} beta diversity calculation")

        return results
    except Exception as e:
        print(f"Error in calculate_beta_diversity: {str(e)}")
        raise


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
    try:
        print("Starting group comparison calculation")
        if not data.metadata or group_by is None:
            print("No metadata or group_by provided, skipping group comparison")
            return None

        groups: Dict[str, List[str]] = {}

        for sample_id in data.sample_ids:
            if sample_id in data.metadata and group_by in data.metadata[sample_id]:
                group_value = data.metadata[sample_id][group_by]

                # Skip NaN values or None values in metadata
                if group_value is None or (
                    isinstance(group_value, float) and np.isnan(group_value)
                ):
                    print(f"Skipping sample {sample_id} due to None/NaN group value")
                    continue

                # Convert group value to string to ensure it works as a dictionary key
                group_value = str(group_value)

                if group_value not in groups:
                    groups[group_value] = []
                groups[group_value].append(sample_id)

        print(f"Found {len(groups)} groups: {list(groups.keys())}")
        if len(groups) < 2:
            print("Less than 2 groups found, skipping group comparison")
            return None

        group_comparison = {}

        for metric, values in diversity_metrics.items():
            print(f"Comparing groups for {metric}")
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

                    print(
                        f"Comparing {group_a} (n={len(group_a_values)}) vs "
                        f"{group_b} (n={len(group_b_values)})"
                    )
                    if len(group_a_values) < 1 or len(group_b_values) < 1:
                        print("Skipping comparison due to insufficient samples")
                        continue

                    # Check for NaN values
                    group_a_values = [v for v in group_a_values if not np.isnan(v)]
                    group_b_values = [v for v in group_b_values if not np.isnan(v)]

                    if len(group_a_values) < 1 or len(group_b_values) < 1:
                        print(
                            "Skipping comparison due to all NaN values after filtering"
                        )
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
            print(f"Completed group comparison for {metric}")

        return group_comparison
    except Exception as e:
        print(f"Error in calculate_group_comparison: {str(e)}")
        raise


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
    try:
        print("Starting diversity analysis")
        print(f"Data: {len(data.sample_ids)} samples, {len(data.feature_ids)} features")
        print(f"Metrics: {metrics}, Group by: {group_by}")

        # Check for NaN values in the counts matrix
        counts_matrix = np.array(data.counts_matrix)
        print(f"Contains NaN values: {np.isnan(counts_matrix).any()}")
        if np.isnan(counts_matrix).any():
            print(f"Number of NaN values: {np.isnan(counts_matrix).sum()}")
            print(f"NaN positions: {np.where(np.isnan(counts_matrix))}")

        alpha_diversity = calculate_alpha_diversity(data, metrics)
        print("Alpha diversity calculation completed")

        beta_diversity = calculate_beta_diversity(data)
        print("Beta diversity calculation completed")

        group_comparison = calculate_group_comparison(data, group_by, alpha_diversity)
        print("Group comparison calculation completed")

        print("Diversity analysis completed successfully")
        return alpha_diversity, beta_diversity, group_comparison
    except Exception as e:
        print(f"Error in run_diversity_analysis: {str(e)}")
        raise
