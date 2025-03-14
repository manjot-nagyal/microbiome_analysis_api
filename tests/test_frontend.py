import json
import os
import sys
from unittest.mock import MagicMock, patch

import pytest

from frontend import (
    perform_correlation_analysis,
    perform_diversity_analysis,
    upload_file,
)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


@pytest.fixture
def sample_data():
    """Fixture providing sample data for tests"""
    return {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["feature1", "feature2"],
        "counts_matrix": [[10, 20], [30, 40]],
    }


@patch("frontend.requests.post")
def test_upload_csv_file(mock_post):
    """Test uploading a CSV file"""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["feature1", "feature2"],
        "counts_matrix": [[10, 20], [30, 40]],
    }
    mock_post.return_value = mock_response

    mock_file = MagicMock()
    mock_file.name = "test.csv"
    mock_file.getvalue.return_value = (
        b"sample_id,feature1,feature2\nsample1,10,20\nsample2,30,40"
    )

    result = upload_file(mock_file, "csv")

    mock_post.assert_called_once()
    assert result == mock_response.json.return_value


@patch("frontend.requests.post")
def test_upload_json_file(mock_post):
    """Test uploading a JSON file"""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["feature1", "feature2"],
        "counts_matrix": [[10, 20], [30, 40]],
    }
    mock_post.return_value = mock_response

    test_data = {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["feature1", "feature2"],
        "counts_matrix": [[10, 20], [30, 40]],
    }

    mock_file = MagicMock()
    mock_file.name = "test.json"
    mock_file.getvalue.return_value = json.dumps(test_data).encode()

    result = upload_file(mock_file, "json")

    mock_post.assert_called_once()
    assert result == mock_response.json.return_value


@patch("frontend.requests.post")
def test_json_transformation(mock_post):
    """Test JSON transformation for alternate formats"""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["taxon1", "taxon2"],
        "counts_matrix": [[0.5, 0.2], [0.1, 0.4]],
    }
    mock_post.return_value = mock_response

    test_data = {
        "sample_ids": ["sample1", "sample2"],
        "taxonomy_ids": ["taxon1", "taxon2"],
        "abundance_matrix": [[0.5, 0.2], [0.1, 0.4]],
    }

    mock_file = MagicMock()
    mock_file.name = "test.json"
    mock_file.getvalue.return_value = json.dumps(test_data).encode()

    result = upload_file(mock_file, "json")

    mock_post.assert_called_once()
    assert result == mock_response.json.return_value


@patch("frontend.requests.post")
def test_perform_diversity_analysis(mock_post):
    """Test the diversity analysis function"""
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = {
        "alpha_diversity_metrics": {"shannon": {"sample1": 0.8, "sample2": 0.9}},
        "beta_diversity_metrics": {"braycurtis": [[0, 0.5], [0.5, 0]]},
        "sample_ids": ["sample1", "sample2"],
    }
    mock_post.return_value = mock_response

    test_data = {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["feature1", "feature2"],
        "counts_matrix": [[10, 20], [30, 40]],
    }
    metrics = ["shannon", "simpson"]

    result = perform_diversity_analysis(test_data, metrics)

    mock_post.assert_called_once()
    assert result == mock_response.json.return_value


@patch("frontend.requests.post")
def test_perform_correlation_analysis(mock_post):
    """Test the correlation analysis function"""
    mock_response1 = MagicMock()
    mock_response1.status_code = 200

    mock_response2 = MagicMock()
    mock_response2.status_code = 200
    mock_response2.json.return_value = {
        "feature_correlation_matrix": [[1.0, 0.8], [0.8, 1.0]],
        "metadata_correlations": {"metadata_age": {"feature1": 0.7}},
    }

    mock_post.side_effect = [mock_response1, mock_response2]

    test_data = {
        "sample_ids": ["sample1", "sample2"],
        "feature_ids": ["feature1", "feature2"],
        "counts_matrix": [[10, 20], [30, 40]],
        "metadata": {"sample1": {"age": 25}, "sample2": {"age": 30}},
    }

    result = perform_correlation_analysis(
        test_data, "spearman", ["metadata_age"], 0.01, 0.2
    )

    assert mock_post.call_count == 2
    assert result == mock_response2.json.return_value
