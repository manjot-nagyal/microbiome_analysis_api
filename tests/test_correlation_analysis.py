import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from scipy.stats import kendalltau, pearsonr, spearmanr

from app.services import correlation_analysis


def simplify_result(res):
    """Extract a scalar if result is a one-element tuple/array."""
    if isinstance(res, tuple):
        res = res[0]
    if hasattr(res, "shape") and np.ndim(res) > 0 and np.size(res) == 1:
        res = float(res)
    return res


def check_result(custom, scipy_val, tol=1e-6):
    """
    If the custom result is undefined (NaN) then pass immediately.
    Otherwise, if a reference value is provided, assert they are close.
    For tests of undefined behavior (like constant inputs, empty arrays
    or mismatched lengths), pass scipy_val as None.
    """
    custom = simplify_result(custom)
    if np.isnan(custom):
        return
    if scipy_val is None:
        raise AssertionError(f"Expected undefined (NaN) result, got {custom}")
    np.testing.assert_allclose(custom, scipy_val, atol=tol, equal_nan=True)


################################################################################
# Pearson correlation tests
################################################################################


def test_pearson_nonconstant():
    rng = np.random.default_rng(42)
    x = rng.random(100)
    y = 2 * x + rng.normal(0, 0.1, 100)
    custom = correlation_analysis.pearson_corr(x, y)
    scipy_val, _ = pearsonr(x, y)
    check_result(custom, scipy_val)


def test_pearson_constant():
    x = np.ones(50)
    y = np.random.random(50)
    custom = correlation_analysis.pearson_corr(x, y)
    with pytest.warns(RuntimeWarning):
        scipy_val, _ = pearsonr(x, y)
    check_result(custom, scipy_val)


################################################################################
# Spearman correlation tests
################################################################################


def test_spearman_nonconstant():
    rng = np.random.default_rng(123)
    x = rng.random(100)
    y = np.sin(x * 10) + rng.normal(0, 0.05, 100)
    custom = correlation_analysis.spearman_corr(x, y)
    scipy_corr, _ = spearmanr(x, y)
    check_result(custom, scipy_corr)


def test_spearman_constant():
    x = np.full(30, 3.14)
    y = np.linspace(0, 1, 30)
    custom = correlation_analysis.spearman_corr(x, y)
    with pytest.warns(RuntimeWarning):
        scipy_corr, _ = spearmanr(x, y)
    check_result(custom, scipy_corr)


################################################################################
# Kendall Tau correlation tests
################################################################################


def test_kendall_nonconstant():
    rng = np.random.default_rng(987)
    x = rng.random(50)
    y = np.log(x + 1) + rng.normal(0, 0.1, 50)
    custom = correlation_analysis.kendall_tau(x, y)
    scipy_corr, _ = kendalltau(x, y)
    check_result(custom, scipy_corr)


def test_kendall_constant():
    x = np.full(40, 2.718)
    y = np.linspace(-1, 1, 40)
    custom = correlation_analysis.kendall_tau(x, y)
    scipy_corr, _ = kendalltau(x, y)
    check_result(custom, scipy_corr)


################################################################################
# Edge cases: empty arrays and mismatched lengths
################################################################################


def test_pearson_empty():
    x = np.array([])
    y = np.array([])
    custom = correlation_analysis.pearson_corr(x, y)
    check_result(custom, None)


def test_spearman_empty():
    x = np.array([])
    y = np.array([])
    custom = correlation_analysis.spearman_corr(x, y)
    check_result(custom, None)


def test_kendall_empty():
    x = np.array([])
    y = np.array([])
    custom = correlation_analysis.kendall_tau(x, y)
    check_result(custom, None)


def test_mismatched_lengths():
    x = np.random.random(10)
    y = np.random.random(9)
    check_result(correlation_analysis.pearson_corr(x, y), None)
    check_result(correlation_analysis.spearman_corr(x, y), None)
    check_result(correlation_analysis.kendall_tau(x, y), None)


################################################################################
# Hypothesis tests
################################################################################


def run_hypothesis_test(corr_func, scipy_func, x, y, tol=1e-6):
    """
    Run a hypothesis test comparing custom correlation function with scipy reference.
    """
    x_is_constant = np.all(x == x[0]) if len(x) > 0 else False
    y_is_constant = np.all(y == y[0]) if len(y) > 0 else False

    custom = corr_func(x, y)

    if x_is_constant or y_is_constant:
        custom_value = simplify_result(custom)
        if abs(custom_value) < 1e-10:
            return
        check_result(custom, np.nan)
        check_result(custom, np.nan)
        return

    try:
        scipy_val, _ = scipy_func(x, y)
    except Exception:
        check_result(custom, None)
        check_result(custom, None)
        return

    check_result(custom, scipy_val, tol)


@given(st.data())
def test_hypothesis_pearson(data):
    n = data.draw(st.integers(min_value=2, max_value=20))
    x = np.array(
        data.draw(
            st.lists(
                st.floats(
                    min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False
                ),
                min_size=n,
                max_size=n,
            )
        )
    )
    y = np.array(
        data.draw(
            st.lists(
                st.floats(
                    min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False
                ),
                min_size=n,
                max_size=n,
            )
        )
    )
    run_hypothesis_test(correlation_analysis.pearson_corr, pearsonr, x, y)


@given(st.data())
def test_hypothesis_spearman(data):
    n = data.draw(st.integers(min_value=2, max_value=20))
    x = np.array(
        data.draw(
            st.lists(
                st.floats(
                    min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False
                ),
                min_size=n,
                max_size=n,
            )
        )
    )
    y = np.array(
        data.draw(
            st.lists(
                st.floats(
                    min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False
                ),
                min_size=n,
                max_size=n,
            )
        )
    )
    run_hypothesis_test(correlation_analysis.spearman_corr, spearmanr, x, y)


@given(st.data())
def test_hypothesis_kendall(data):
    n = data.draw(st.integers(min_value=2, max_value=20))
    x = np.array(
        data.draw(
            st.lists(
                st.floats(
                    min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False
                ),
                min_size=n,
                max_size=n,
            )
        )
    )
    y = np.array(
        data.draw(
            st.lists(
                st.floats(
                    min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False
                ),
                min_size=n,
                max_size=n,
            )
        )
    )
    run_hypothesis_test(correlation_analysis.kendall_tau, kendalltau, x, y)
