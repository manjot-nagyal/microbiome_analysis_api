import pytest
from fastapi.testclient import TestClient

from app.main import app  # Import your FastAPI app
from app.models.schemas import DiversityAnalysisRequest, MicrobiomeData

client = TestClient(app)


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


def test_upload_no_file():
    response = client.post("/upload/")
    assert response.status_code == 400
    assert response.json() == {"detail": "No file uploaded"}


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
    response = client.post(
        "/upload/", files={"file": ("test.csv", csv_content, "text/csv")}
    )
    assert response.status_code == 400
    assert "Error parsing CSV file" in response.json()["detail"]


def test_upload_csv_with_parsing_error():
    csv_content = "sample_id,feature1\nsample1,10\nsample2,not_a_number"
    response = client.post(
        "/upload/", files={"file": ("test.csv", csv_content, "text/csv")}
    )
    assert response.status_code == 400
    assert "Error parsing CSV file" in response.json()["detail"]


def test_analyze_diversity(sample_data):
    # Sample request data
    request = DiversityAnalysisRequest(
        metrics=["shannon", "simpson", "pielou", "chao1"]
    )

    # Sample microbiome data
    data = sample_data

    payload = {"request": {"metrics": request.metrics}, "data": data.model_dump()}

    # Make the API request
    response = client.post("/diversity/analyze", json=payload)

    print(response.json())

    # Assert the response
    assert response.status_code == 200

    # Assert response structure
    assert "alpha_diversity_metrics" in response.json()
    assert "beta_diversity_metrics" in response.json()
    assert "group_comparison_metrics" in response.json()
    assert "sample_ids" in response.json()

    # Assert response values
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
