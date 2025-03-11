import numpy as np
import pytest
from scipy.spatial import distance


from app.models.schemas import MicrobiomeData
from app.services.diversity_analysis import (
    calculate_alpha_diversity,
    calculate_beta_diversity,
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
            "sample5": {"group": "C", "age": 25, "sex": "M"},
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

    for metric, matrix in results.items():
        # Check that there are results for every metric
        assert metric in results
        # Account for the size of the matrix matrches that of the number of features
        square_matrix = distance.squareform(results[metric])
        assert len(square_matrix) == len(sample_data.feature_ids)
        for row in square_matrix:
            assert len(row) == len(sample_data.feature_ids)

        # Check that the matrix is symmetric
        assert np.allclose(square_matrix, square_matrix.T)

        # Check that the diagonal is 0 (distance to self)
        assert np.allclose(
            np.diag(square_matrix), np.zeros(len(sample_data.feature_ids))
        )
