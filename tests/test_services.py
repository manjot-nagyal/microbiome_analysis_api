import numpy as np
import pytest

from app.models.schemas import MicrobiomeData
from app.services.diversity_analysis import (
    calculate_alpha_diversity,
    calculate_beta_diversity,
    calculate_group_comparison,
)


# Sample test data
@pytest.fixture
def sample_data():
    return MicrobiomeData(
        sample_ids=["sample1", "sample2", "sample3", "sample4", "sample5"],
        feature_ids=["taxon1", "taxon2", "taxon3", "taxon4"],
        counts_matrix=[
            [0.5, 0.2, 0.1, 0.1, 0.1],
            [0.1, 0.4, 0.2, 0.2, 0.1],
            [0.3, 0.3, 0.2, 0.1, 0.1],
            [0.2, 0.2, 0.3, 0.2, 0.1],
        ],
        metadata={
            "sample1": {"group": "A", "age": 25, "sex": "M"},
            "sample2": {"group": "A", "age": 30, "sex": "F"},
            "sample3": {"group": "B", "age": 28, "sex": "M"},
            "sample4": {"group": "B", "age": 32, "sex": "F"},
            "sample5": {"group": "A", "age": 25, "sex": "M"},
        },
    )


def test_calculate_alpha_diversity(sample_data):
    """
    Test the alpha diversity calculation
    """

    metrics = ["shannon", "simpson", "pielou", "chao1"]
    results = calculate_alpha_diversity(sample_data, metrics)

    # Check that there are results for every sample
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
        # Check that there are results for every metric
        assert metric in results

        # Account for the size of the matrix matrches that of the number of samples
        for row in square_matrix:
            assert len(row) == len(sample_data.sample_ids)

        square_matrix_np = np.array(square_matrix)
        # Check that the matrix is symmetric
        assert np.allclose(square_matrix_np, square_matrix_np.T)

        # Check that the diagonal is 0 (distance to self)
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
