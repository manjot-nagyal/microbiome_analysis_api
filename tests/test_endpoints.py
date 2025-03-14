import json

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.models.schemas import DiversityAnalysisRequest, MicrobiomeData

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

    # Test invalid diversity metric
    client.post(
        "/upload/",
        data={"json_data": json.dumps(json_data), "normalize": "true"},
    )

    response = client.post(
        "/diversity/analyze",
        json={"request": {"metrics": ["invalid_metric"]}, "data": json_data},
    )
    assert response.status_code == 400
