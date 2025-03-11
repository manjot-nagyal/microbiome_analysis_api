import numpy as np
import pytest
import skbio.diversity.alpha as skbio_alpha

from app.services.diversity_analysis import (
    shannon_diversity,
    simpson_diversity,
    pielou_evenness,
    chao1_estimator,
)


# Helper functions:
def assert_alpha_metric(result, expected, tolerance=0.001):
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
        """Check for uniform distribution"""
        uniform = np.ones(10)
        expected = np.log(10)

        result = shannon_diversity(uniform)
        assert_alpha_metric(result, expected)

    def test_shannon_diversity_single_value(self):
        """Check for single value"""
        single_value = np.array([1])
        expected = np.nan

        result = shannon_diversity(single_value)
        assert_alpha_metric(result, expected)

    def test_shannon_diversity_zero_values(self):
        """Check for zero values"""
        zero_values = np.array([0, 0, 0])
        expected = np.nan

        result = shannon_diversity(zero_values)
        assert_alpha_metric(result, expected)

    def test_shannon_diversity_negative_values(self):
        """Check for negative values"""
        negative_values = np.array([-1, -2, -3])
        with pytest.raises(ValueError):
            shannon_diversity(negative_values)

    def test_shannon_diversity_large_values(self):
        """Check for large values"""
        large_values = np.array([1e6, 1e6, 1e6])
        expected = skbio_alpha.shannon(large_values)

        result = shannon_diversity(large_values)
        assert_alpha_metric(result, expected)

    def test_shannon_diversity_mixed_values(self):
        """Check for mixed values"""
        mixed_values = np.array([1, 2, 3, 4, 5])
        expected = skbio_alpha.shannon(mixed_values)

        result = shannon_diversity(mixed_values)
        assert_alpha_metric(result, expected)

    ###########################################################################
    # Simpson diversity tests
    ###########################################################################
    def test_simpson_diversity_uniform(self):
        """Check for uniform distribution"""
        uniform = np.ones(10)
        expected = skbio_alpha.simpson(uniform)

        result = simpson_diversity(uniform)
        assert_alpha_metric(result, expected)

    def test_simpson_diversity_single_value(self):
        """Check for single value"""
        single_value = np.array([1])
        expected = np.nan

        result = simpson_diversity(single_value)
        assert_alpha_metric(result, expected)

    def test_simpson_diversity_zero_values(self):
        """Check for zero values"""
        zero_values = np.array([0, 0, 0])
        expected = np.nan

        result = simpson_diversity(zero_values)
        assert_alpha_metric(result, expected)

    def test_simpson_diversity_negative_values(self):
        """Check for negative values"""
        negative_values = np.array([-1, -2, -3])
        with pytest.raises(ValueError):
            simpson_diversity(negative_values)

    def test_simpson_diversity_mixed_values(self):
        """Check for mixed values"""
        mixed_values = np.array([1, 2, 3, 4, 5])
        expected = skbio_alpha.simpson(mixed_values)

        result = simpson_diversity(mixed_values)
        assert_alpha_metric(result, expected)

    ###########################################################################
    # Pielou's evenness tests
    ###########################################################################
    def test_pielou_evenness_uniform(self):
        """Check for uniform distribution"""
        uniform = np.ones(10)
        expected = skbio_alpha.pielou_e(uniform)

        result = pielou_evenness(uniform)
        assert_alpha_metric(result, expected)

    def test_pielou_evenness_single_value(self):
        """Check for single value"""
        single_value = np.array([1, 0, 0, 0, 0])
        expected = skbio_alpha.pielou_e(single_value)

        result = pielou_evenness(single_value)
        assert_alpha_metric(result, expected)

    def test_pielou_evenness_zero_values(self):
        """Check for zero values"""
        zero_values = np.array([0, 0, 0])
        expected = skbio_alpha.pielou_e(zero_values)

        result = pielou_evenness(zero_values)
        assert_alpha_metric(result, expected)

    def test_pielou_evenness_negative_values(self):
        """Check for negative values"""
        negative_values = np.array([-1, -2, -3])
        with pytest.raises(ValueError):
            pielou_evenness(negative_values)

    def test_pielou_evenness_mixed_values(self):
        """Check for mixed values"""
        mixed_values = np.array([1, 2, 3, 4, 5])
        expected = skbio_alpha.pielou_e(mixed_values)

        result = pielou_evenness(mixed_values)
        assert_alpha_metric(result, expected)

    ###########################################################################
    # Chao1 estimator tests
    ###########################################################################

    def test_chao1_estimator_uniform(self):
        """Check for uniform distribution"""
        uniform = np.ones(10)
        expected = skbio_alpha.chao1(uniform)

        result = chao1_estimator(uniform)
        assert_alpha_metric(result, expected)

    def test_chao1_estimator_single_value(self):
        """Check for single value"""
        single_value = np.array([1])
        expected = skbio_alpha.chao1(single_value)

        result = chao1_estimator(single_value)
        assert_alpha_metric(result, expected)

    def test_chao1_estimator_doubletons(self):
        """Check for doubletons"""
        doubletons = np.array([2, 2, 2, 2, 2])
        expected = skbio_alpha.chao1(doubletons)

        result = chao1_estimator(doubletons)
        assert_alpha_metric(result, expected)

    def test_chao1_estimator_zero_values(self):
        """Check for zero values"""
        zero_values = np.array([0, 0, 0])
        expected = np.nan

        result = chao1_estimator(zero_values)
        assert_alpha_metric(result, expected)

    def test_chao1_estimator_negative_values(self):
        """Check for negative values"""
        negative_values = np.array([-1, -2, -3])
        with pytest.raises(ValueError):
            chao1_estimator(negative_values)

    def test_chao1_estimator_mixed_values(self):
        """Check for mixed values"""
        mixed_values = np.array([1, 2, 3, 4, 5])
        expected = skbio_alpha.chao1(mixed_values)

        result = chao1_estimator(mixed_values)
        assert_alpha_metric(result, expected)

    def test_edge_cases(self):
        # Empty array
        empty = np.array([])
        assert np.isnan(shannon_diversity(empty))
        assert np.isnan(simpson_diversity(empty))
        assert np.isnan(pielou_evenness(empty))
        assert np.isnan(chao1_estimator(empty))
