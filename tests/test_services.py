import numpy as np
import pytest

from app.models.schemas import MicrobiomeData
from app.services.correlation_analysis import (
    calculate_taxon_correlations,
    filter_taxa,
    run_correlation_analysis,
)
from app.services.diversity_analysis import (
    calculate_alpha_diversity,
    calculate_beta_diversity,
    calculate_group_comparison,
)
from app.utils.parsers import normalize_counts


@pytest.fixture
def sample_data():
    return MicrobiomeData(
        sample_ids=["sample1", "sample2", "sample3", "sample4", "sample5"],
        feature_ids=["taxon1", "taxon2", "taxon3", "taxon4"],
        counts_matrix=[
            [0.5, 0.2, 0.1, 0.1],
            [0.1, 0.4, 0.2, 0.2],
            [0.3, 0.3, 0.2, 0.1],
            [0.2, 0.2, 0.3, 0.2],
            [0.1, 0.1, 0.1, 0.1],
        ],
        metadata={
            "sample1": {"group": "A", "age": 25, "sex": "M"},
            "sample2": {"group": "A", "age": 30, "sex": "F"},
            "sample3": {"group": "B", "age": 28, "sex": "M"},
            "sample4": {"group": "B", "age": 32, "sex": "F"},
            "sample5": {"group": "A", "age": 25, "sex": "M"},
        },
    )


@pytest.fixture
def sample_data_with_nan():
    """Create sample data with NaN and empty values"""
    return MicrobiomeData(
        sample_ids=["sample1", "sample2", "sample3", "sample4", "sample5"],
        feature_ids=["taxon1", "taxon2", "taxon3", "taxon4"],
        counts_matrix=[
            [0.5, np.nan, 0.1, 0.1],  # NaN in count matrix
            [0.1, 0.4, 0.2, 0.2],
            [0.3, 0.3, np.nan, 0.1],  # NaN in count matrix
            [0.2, 0.2, 0.3, 0.2],
            [0.1, 0.1, 0.1, 0.1],
        ],
        metadata={
            "sample1": {"group": "A", "age": 25, "sex": "M"},
            "sample2": {"group": "A", "age": np.nan, "sex": "F"},  # NaN in metadata
            "sample3": {"group": "B", "age": 28, "sex": "M"},
            "sample4": {"group": "B", "age": 32, "sex": ""},  # Empty string in metadata
            "sample5": {"group": "C", "age": 25, "sex": "M"},
        },
    )


def test_calculate_alpha_diversity(sample_data):
    """
    Test the alpha diversity calculation
    """

    metrics = ["shannon", "simpson", "pielou", "chao1"]
    results = calculate_alpha_diversity(sample_data, metrics)

    for _, values in results.items():
        assert isinstance(values, dict)
        assert len(values) == len(sample_data.sample_ids)
        for sample_id in sample_data.sample_ids:
            assert sample_id in values


def test_calculate_beta_diversity(sample_data):
    """
    Test the beta diversity calculation
    """
    metrics = ["braycurtis", "jaccard"]
    results = calculate_beta_diversity(sample_data, metrics)

    for metric, square_matrix in results.items():

        assert metric in results

        for row in square_matrix:
            assert len(row) == len(sample_data.sample_ids)

        square_matrix_np = np.array(square_matrix)

        assert np.allclose(square_matrix_np, square_matrix_np.T)

        assert np.allclose(
            np.diag(square_matrix), np.zeros(len(sample_data.sample_ids))
        )


def test_calculate_group_comparison(sample_data):
    """
    Test the group comparison calculation
    """
    alpha_diversity_metrics = calculate_alpha_diversity(sample_data)
    results = calculate_group_comparison(sample_data, "group", alpha_diversity_metrics)

    assert "shannon" in results
    assert "simpson" in results
    assert "pielou" in results
    assert "chao1" in results

    for comparison in results.values():
        assert "A_vs_B" in comparison

        comparison = comparison["A_vs_B"]

        assert "p_value" in comparison
        assert "statistic" in comparison
        assert "significant" in comparison
        assert "A_mean" in comparison
        assert "B_mean" in comparison


def test_normalize_counts(sample_data):
    """Test that abundance data is properly normalized."""
    non_normalized = MicrobiomeData(
        sample_ids=["sample1", "sample2"],
        feature_ids=["taxon1", "taxon2", "taxon3"],
        counts_matrix=[[10, 20, 70], [50, 30, 20]],
    )

    normalized = normalize_counts(non_normalized)

    # Check that each row sums to 1.0
    for row in normalized.counts_matrix:
        assert abs(sum(row) - 1.0) < 1e-10

    # Check that relative proportions are maintained
    assert normalized.counts_matrix[0][0] == 0.1
    assert normalized.counts_matrix[0][1] == 0.2
    assert normalized.counts_matrix[0][2] == 0.7
    assert normalized.counts_matrix[1][0] == 0.5
    assert normalized.counts_matrix[1][1] == 0.3
    assert normalized.counts_matrix[1][2] == 0.2


def test_filter_taxa(sample_data):
    """Test filtering taxa based on abundance and prevalence."""
    test_data = MicrobiomeData(
        sample_ids=["sample1", "sample2", "sample3", "sample4"],
        feature_ids=["common_taxon", "rare_taxon", "low_prev_taxon"],
        counts_matrix=[
            [0.5, 0.001, 0.1],
            [0.6, 0.002, 0.0],
            [0.7, 0.003, 0.0],
            [0.8, 0.0, 0.0],
        ],
    )

    filtered_matrix, filtered_taxa = filter_taxa(
        test_data, min_abundance=0.01, min_prevalence=0.5
    )

    assert len(filtered_taxa) == 1
    assert filtered_taxa[0] == "common_taxon"
    assert filtered_matrix.shape == (4, 1)

    filtered_matrix, filtered_taxa = filter_taxa(
        test_data, min_abundance=0.01, min_prevalence=0.2
    )
    assert len(filtered_taxa) == 2
    assert "common_taxon" in filtered_taxa
    assert "low_prev_taxon" in filtered_taxa

    filtered_matrix, filtered_taxa = filter_taxa(
        test_data, min_abundance=0.001, min_prevalence=0.5
    )
    assert len(filtered_taxa) == 2
    assert "common_taxon" in filtered_taxa
    assert "rare_taxon" in filtered_taxa


def test_taxon_correlations(sample_data):
    """Test calculating correlations between taxa."""
    filtered_matrix, filtered_taxa = filter_taxa(sample_data)

    corr_dict, p_dict = calculate_taxon_correlations(
        filtered_matrix, filtered_taxa, method="spearman"
    )

    for taxon1 in filtered_taxa:
        assert taxon1 in corr_dict
        assert taxon1 in p_dict

        for taxon2 in filtered_taxa:
            if taxon1 != taxon2:
                assert taxon2 in corr_dict[taxon1]
                assert taxon2 in p_dict[taxon1]

                # Correlation value should be between -1 and 1
                assert -1.0 <= corr_dict[taxon1][taxon2] <= 1.0

                # P-value should be between 0 and 1
                assert 0.0 <= p_dict[taxon1][taxon2] <= 1.0


def test_run_correlation_analysis(sample_data):
    """Test the full correlation analysis workflow."""
    taxon_corr, meta_corr, p_values, filtered_taxa = run_correlation_analysis(
        sample_data,
        correlation_method="spearman",
        metadata_columns=["age", "group"],
        min_abundance=0.01,
        min_prevalence=0.1,
    )

    assert taxon_corr
    assert meta_corr
    assert p_values
    assert filtered_taxa

    for taxon in filtered_taxa:
        assert taxon in meta_corr
        assert "age" in meta_corr[taxon]
