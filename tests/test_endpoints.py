import json

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.models.schemas import (
    CorrelationAnalysisRequest,
    DiversityAnalysisRequest,
    MicrobiomeData,
)

client = TestClient(app)


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
            "sample5": {"group": "C", "age": 25, "sex": "M"},
        },
    )


@pytest.fixture
def json_data():
    return {
        "sample_ids": ["sample1", "sample2", "sample3", "sample4", "sample5"],
        "feature_ids": ["taxon1", "taxon2", "taxon3", "taxon4"],
        "counts_matrix": [
            [0.5, 0.2, 0.1, 0.1],
            [0.1, 0.4, 0.2, 0.2],
            [0.3, 0.3, 0.2, 0.1],
            [0.2, 0.2, 0.3, 0.2],
            [0.1, 0.1, 0.1, 0.1],
        ],
        "metadata": {
            "sample1": {"group": "A", "age": 25, "sex": "M"},
            "sample2": {"group": "A", "age": 30, "sex": "F"},
            "sample3": {"group": "B", "age": 28, "sex": "M"},
            "sample4": {"group": "B", "age": 32, "sex": "F"},
            "sample5": {"group": "C", "age": 25, "sex": "M"},
        },
    }


@pytest.fixture
def sample_data_with_nan():
    return MicrobiomeData(
        sample_ids=["sample1", "sample2", "sample3", "sample4", "sample5"],
        feature_ids=["taxon1", "taxon2", "taxon3", "taxon4"],
        counts_matrix=[
            [0.5, 0.2, 0.0, 0.1],
            [0.1, 0.4, 0.0, 0.2],
            [0.3, 0.3, 0.2, 0.0],
            [0.2, 0.0, 0.3, 0.2],
            [0.1, 0.1, 0.0, 0.1],
        ],
        metadata={
            "sample1": {"group": "A", "age": 25, "sex": "M"},
            "sample2": {"group": "A", "age": 30, "sex": "F"},
            "sample3": {"group": "B", "age": 28, "sex": "M"},
            "sample4": {"group": "B", "age": 32, "sex": "F"},
            "sample5": {"group": "C", "age": 25, "sex": "M"},
        },
    )


def test_upload_no_file():
    response = client.post("/upload/")
    assert response.status_code == 400
    assert response.json() == {"detail": "No file uploaded or JSON data provided"}


def test_upload_invalid_file_type():
    response = client.post(
        "/upload/",
        files={
            "file": (
                "test.txt",
                "This is not a CSV file",
                "text/plain",
            )
        },
    )
    assert response.status_code == 400
    assert "Unsupported file type" in response.json()["detail"]


def test_upload_valid_csv():
    csv_content = "sample_id,feature1,feature2\nsample1,10,20\nsample2,30,40"

    from io import BytesIO

    file_like = BytesIO(csv_content.encode("utf-8"))

    response = client.post(
        "/upload/", files={"file": ("test.csv", file_like, "text/csv")}
    )
    response_data = response.json()
    assert response.status_code == 200
    assert "sample_ids" in response_data
    assert response_data["sample_ids"] == ["sample1", "sample2"]


def test_upload_csv_with_parsing_error():
    csv_content = "sample_id,feature1\nsample1,10\nsample2,not_a_number"
    response = client.post(
        "/upload/", files={"file": ("test.csv", csv_content, "text/csv")}
    )
    assert response.status_code == 400
    assert "Error parsing CSV file" in response.json()["detail"]


def test_upload_json_data_direct(json_data):
    """Test uploading microbiome data as JSON."""
    response = client.post(
        "/upload/",
        data={"json_data": json.dumps(json_data), "normalize": "true"},
    )
    assert response.status_code == 200

    data = response.json()
    assert data["sample_ids"] == json_data["sample_ids"]
    assert data["feature_ids"] == json_data["feature_ids"]
    assert len(data["counts_matrix"]) == len(json_data["counts_matrix"])


def test_upload_json_data_file(json_data):
    """Test uploading microbiome data as JSON."""
    response = client.post(
        "/upload/",
        files={"file": ("test.json", json.dumps(json_data), "application/json")},
    )
    assert response.status_code == 200
    assert "sample_ids" in response.json()
    assert response.json()["sample_ids"] == [
        "sample1",
        "sample2",
        "sample3",
        "sample4",
        "sample5",
    ]


def test_analyze_diversity(sample_data):
    # Sample request data
    request = DiversityAnalysisRequest(
        metrics=["shannon", "simpson", "pielou", "chao1"]
    )
    data = sample_data
    payload = {"request": {"metrics": request.metrics}, "data": data.model_dump()}

    response = client.post("/diversity/analyze", json=payload)

    assert response.status_code == 200

    assert "alpha_diversity_metrics" in response.json()
    assert "beta_diversity_metrics" in response.json()
    assert "group_comparison_metrics" in response.json()
    assert "sample_ids" in response.json()

    assert len(response.json()["alpha_diversity_metrics"]["shannon"]) == len(
        data.sample_ids
    )
    assert len(response.json()["alpha_diversity_metrics"]["simpson"]) == len(
        data.sample_ids
    )
    assert len(response.json()["alpha_diversity_metrics"]["pielou"]) == len(
        data.sample_ids
    )
    assert len(response.json()["alpha_diversity_metrics"]["chao1"]) == len(
        data.sample_ids
    )
    assert len(response.json()["beta_diversity_metrics"]["braycurtis"]) == len(
        data.sample_ids
    )
    assert len(response.json()["beta_diversity_metrics"]["jaccard"]) == len(
        data.sample_ids
    )
    assert len(response.json()["sample_ids"]) == len(data.sample_ids)


def test_error_handling(json_data):
    """Test error handling in the API."""

    client.post(
        "/upload/",
        data={"json_data": json.dumps(json_data), "normalize": "true"},
    )

    response = client.post(
        "/diversity/analyze",
        json={"request": {"metrics": ["invalid_metric"]}, "data": json_data},
    )
    assert response.status_code == 400


def test_correlation_analyze_basic(sample_data):
    """Test basic correlation analysis with default parameters."""

    upload_response = client.post(
        "/upload/",
        data={"json_data": json.dumps(sample_data.model_dump()), "normalize": "true"},
    )
    assert upload_response.status_code == 200

    request = CorrelationAnalysisRequest()

    response = client.post(
        "/correlation/analyze", json={"request": request.model_dump()}
    )

    assert response.status_code == 200

    result = response.json()
    assert "taxon_correlations" in result
    assert "p_values" in result
    assert "filtered_taxa" in result

    for taxon in sample_data.feature_ids:
        if taxon in result["filtered_taxa"]:
            assert taxon in result["taxon_correlations"]
            assert len(result["taxon_correlations"][taxon]) > 0


def test_correlation_analyze_methods(sample_data):
    """Test correlation analysis with different correlation methods."""

    upload_response = client.post(
        "/upload/",
        data={"json_data": json.dumps(sample_data.model_dump()), "normalize": "true"},
    )
    assert upload_response.status_code == 200

    for method in ["pearson", "spearman", "kendall"]:
        request = CorrelationAnalysisRequest(correlation_method=method)

        response = client.post(
            "/correlation/analyze", json={"request": request.model_dump()}
        )

        assert response.status_code == 200
        result = response.json()

        assert "taxon_correlations" in result
        assert "p_values" in result
        assert "filtered_taxa" in result


def test_correlation_analyze_with_metadata(sample_data):
    """Test correlation analysis including metadata variables."""

    upload_response = client.post(
        "/upload/",
        data={"json_data": json.dumps(sample_data.model_dump()), "normalize": "true"},
    )
    assert upload_response.status_code == 200

    request = CorrelationAnalysisRequest(
        correlation_method="spearman", metadata_columns=["age", "sex"]
    )

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 200
    result = response.json()

    assert "metadata_correlations" in result

    for taxon in result["filtered_taxa"]:
        assert (
            taxon in result["metadata_correlations"]
        ), f"Taxon {taxon} should have metadata correlations"
        metadata_corrs = result["metadata_correlations"][taxon]
        assert any(col in metadata_corrs for col in ["age", "sex"])


def test_correlation_analyze_filtering(sample_data):
    """Test correlation analysis with different abundance and prevalence thresholds."""

    upload_response = client.post(
        "/upload/",
        data={"json_data": json.dumps(sample_data.model_dump()), "normalize": "true"},
    )
    assert upload_response.status_code == 200

    request = CorrelationAnalysisRequest(min_abundance=0.3, min_prevalence=0.6)

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 200
    result = response.json()

    assert len(result["filtered_taxa"]) <= len(sample_data.feature_ids)

    request = CorrelationAnalysisRequest(min_abundance=0.01, min_prevalence=0.01)

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 200
    result = response.json()

    assert len(result["filtered_taxa"]) >= 1


def test_correlation_analyze_with_zeros_and_nans(sample_data_with_nan):
    """Test correlation analysis with data containing zeros and potential NaN values."""

    upload_response = client.post(
        "/upload/",
        data={
            "json_data": json.dumps(sample_data_with_nan.model_dump()),
            "normalize": "true",
        },
    )
    assert upload_response.status_code == 200

    request = CorrelationAnalysisRequest(
        correlation_method="spearman", min_abundance=0.01, min_prevalence=0.01
    )

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 200

    result = response.json()
    assert "taxon_correlations" in result
    assert "p_values" in result
    assert "filtered_taxa" in result


def test_correlation_analyze_invalid_method():
    """Test correlation analysis with an invalid correlation method."""

    data = MicrobiomeData(
        sample_ids=["sample1", "sample2"],
        feature_ids=["taxon1", "taxon2"],
        counts_matrix=[[0.5, 0.5], [0.5, 0.5]],
        metadata={"sample1": {"group": "A"}, "sample2": {"group": "B"}},
    )

    upload_response = client.post(
        "/upload/",
        data={"json_data": json.dumps(data.model_dump()), "normalize": "true"},
    )
    assert upload_response.status_code == 200

    request = CorrelationAnalysisRequest(correlation_method="invalid_method")

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 400
    assert "Unsupported correlation method" in response.json()["detail"]


def test_correlation_analyze_invalid_metadata_columns(sample_data):
    """Test correlation analysis with invalid metadata columns."""

    upload_response = client.post(
        "/upload/",
        data={"json_data": json.dumps(sample_data.model_dump()), "normalize": "true"},
    )
    assert upload_response.status_code == 200

    request = CorrelationAnalysisRequest(metadata_columns=["nonexistent_column"])

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 400
    assert "Unknown metadata columns" in response.json()["detail"]


def test_correlation_analyze_no_data():
    """Test correlation analysis without first uploading data."""

    client.post("/upload/current")

    request = CorrelationAnalysisRequest()

    response = client.post("/correlation/analyze", json=request.model_dump())

    assert response.status_code == 400
    assert "No microbiome data available" in response.json()["detail"]
