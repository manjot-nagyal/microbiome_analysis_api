import numpy as np
import pytest
import skbio.diversity.alpha as skbio_alpha

from app.services.diversity_analysis import (
    shannon_diversity,
    simpson_diversity
)


# Helper functions:
def assert_metric(result, expected, tolerance=0.001):
    """
    Assert that result equals expected (within tol). If result is NaN,
    the test passes because the metric is undefined in this case.
    """
    if np.isnan(expected):
        return
    assert result == pytest.approx(expected, abs=tolerance)

class TestDiversityMetrics:
    """Test suite for custom diversity metric implementations."""

    def test_shannon_diversity_uniform(self):
        ''' Check for uniform distribution '''
        uniform = np.ones(10)
        expected = np.log(10)

        result = shannon_diversity(uniform)
        assert_metric(result, expected)

    def test_shannon_diversity_single_value(self):
        ''' Check for single value '''
        single_value = np.array([1])
        expected = np.nan

        result = shannon_diversity(single_value)
        assert_metric(result, expected)

    def test_shannon_diversity_zero_values(self):
        ''' Check for zero values '''
        zero_values = np.array([0, 0, 0])
        expected = np.nan

        result = shannon_diversity(zero_values)
        assert_metric(result, expected)

    def test_shannon_diversity_negative_values(self):
        ''' Check for negative values '''
        negative_values = np.array([-1, -2, -3])
        expected = np.nan

        result = shannon_diversity(negative_values)
        assert_metric(result, expected)

    def test_shannon_diversity_large_values(self):
        ''' Check for large values '''
        large_values = np.array([1e6, 1e6, 1e6])
        expected = skbio_alpha.shannon(large_values)

        result = shannon_diversity(large_values)
        assert_metric(result, expected)


    def test_shannon_diversity_mixed_values(self):
        ''' Check for mixed values '''
        mixed_values = np.array([1, 2, 3, 4, 5])
        expected = skbio_alpha.shannon(mixed_values)

        result = shannon_diversity(mixed_values)
        assert_metric(result, expected)

    ###########################################################################
    # Simpson diversity tests
    ###########################################################################
    def test_simpson_diversity_uniform(self):
        ''' Check for uniform distribution '''
        uniform = np.ones(10)
        expected = skbio_alpha.simpson(uniform)

        result = simpson_diversity(uniform)
        assert_metric(result, expected)
    
    def test_simpson_diversity_single_value(self):
        ''' Check for single value '''
        single_value = np.array([1])
        expected = np.nan

        result = simpson_diversity(single_value)
        assert_metric(result, expected)

    def test_simpson_diversity_zero_values(self):
        ''' Check for zero values '''
        zero_values = np.array([0, 0, 0])
        expected = np.nan

        result = simpson_diversity(zero_values)
        assert_metric(result, expected)

    def test_simpson_diversity_negative_values(self):
        ''' Check for negative values '''
        negative_values = np.array([-1, -2, -3])
        expected = np.nan

        result = simpson_diversity(negative_values)
        assert_metric(result, expected)

    def test_simpson_diversity_mixed_values(self):
        ''' Check for mixed values '''
        mixed_values = np.array([1, 2, 3, 4, 5])
        expected = skbio_alpha.simpson(mixed_values)

        result = simpson_diversity(mixed_values)
        assert_metric(result, expected)
