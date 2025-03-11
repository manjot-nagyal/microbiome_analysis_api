from fastapi.testclient import TestClient
from app.main import app  # Import your FastAPI app

client = TestClient(app)


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
