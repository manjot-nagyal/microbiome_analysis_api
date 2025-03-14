# Microbiome Data Analysis API

A containerized REST API for analyzing microbiome data, built with FastAPI. This API allows users to upload microbiome abundance data in JSON or CSV format and perform various diversity and correlation analyses. 


## Features

- **Data Upload**: Support for CSV and JSON file formats
- **Diversity Analysis**: Calculate alpha diversity metrics (Shannon, Simpson, Pielou, Chao1) and beta diversity metrics
- **Correlation Analysis**: Analyze correlations between taxa and metadata variables
- **Interactive Visualizations**: Visualize analysis results with interactive plots

## Installation

### Prerequisites

- Docker and Docker Compose (recommended)
- Python 3.12+ (for local development)
- uv package manager (for local development)

### Option 1: Using Docker (Recommended)

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd microbiome_analysis_api
   ```

2. Build and run the containers:
   ```bash
   docker-compose up -d
   ```

3. Access the application:
   - Frontend: http://localhost:8501
   - Backend API: http://localhost:8000
   - API Documentation (Swagger UI): http://localhost:8000/docs
   - API Documentation (ReDoc): http://localhost:8000/redoc

### Option 2: Local Development

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd microbiome_analysis_api
   ```

2. Install dependencies using uv:
   ```bash
   uv pip install -e ".[dev]"
   ```

3. Run the backend:
   ```bash
   uv run uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
   ```

4. Run the frontend (in a separate terminal):
   ```bash
   uv run streamlit run frontend.py -- --server.port=8501 --server.address=0.0.0.0
   ```

5. Access the application:
   - Frontend: http://localhost:8501
   - Backend API: http://localhost:8000
   - API Documentation (Swagger UI): http://localhost:8000/docs
   - API Documentation (ReDoc): http://localhost:8000/redoc

## Usage

### Data Format Requirements

#### CSV Format
- First column: sample_id
- Remaining Columns: feature_ids (e.g. taxa names)
- Remaining cells: feature counts
- Optional: metadata columns with the 'metadata_' prefix

Example:
| sample_id | taxon1  | taxon2  | metadata\_group | metadata\_treatment |
| --------  | ------- | ------- | --------------- | ------------------- |
| sample1   | 5.2     | 7.8     | A               | control             |
| sample2   | 10.1    | 12.3    | B               | treatment           |


#### JSON Format
```json
{
  "sample_ids": ["sample1", "sample2"],
  "feature_ids": ["feature1", "feature2", "feature3"],
  "counts_matrix": [
    [10, 20, 30],
    [15, 25, 35]
  ],
  "metadata": {
    "sample1": {"group": "A", "ph": 7.2},
    "sample2": {"group": "B", "ph": 6.8}
  }
}
```

Sample data of both formats (JSON and CSV) is available in the `sample_data` directory. 

### Using the Dashboard

1. **Upload Data**:
   - Select the file format (CSV or JSON)
   - Upload your microbiome data file
   - Preview the uploaded data

2. **Diversity Analysis**:
   - Select the diversity metrics to calculate
   - View interactive visualizations of alpha and beta diversity

3. **Correlation Analysis**:
   - Set correlation parameters (method, minimum abundance, minimum prevalence)
   - Select metadata columns for correlation (if available)
   - View correlation heatmaps and network visualizations

## API Endpoints

### API Documentation
- Interactive API docs (Swagger UI): http://localhost:8000/docs
- Alternative API docs (ReDoc): http://localhost:8000/redoc

### Upload Endpoints
- `POST /upload/`: Upload microbiome data (CSV or JSON)
- `GET /upload/current`: Get the currently loaded microbiome data
- `POST /upload/current`: Reset the currently loaded microbiome data

### Diversity Analysis Endpoints
- `POST /diversity/analyze`: Analyze diversity of microbiome data

### Correlation Analysis Endpoints
- `POST /correlation/analyze`: Analyze correlations in microbiome data

## Development

### Running Tests

Using Docker:
```bash
docker build --target test .
```

Using uv:
```bash
uv run pytest -v
```

### Project Structure

```
microbiome_analysis_api/
├── app/                      # Backend API
│   ├── api/                  # API endpoints
│   │   └── endpoints/        # Endpoint modules
│   ├── models/               # Data models
│   ├── services/             # Analysis services
│   └── utils/                # Utility functions
├── tests/                    # Test suite
├── frontend.py               # Streamlit frontend
├── Dockerfile                # Multi-stage Docker build
├── docker-compose.yml        # Docker Compose configuration
└── pyproject.toml            # Project dependencies
```

## Technical Details

### Backend (FastAPI)
- Built with FastAPI for high-performance API endpoints
- Uses Pydantic for data validation
- Implements microbiome analysis algorithms for diversity and correlation analysis

### Frontend (Streamlit)
- Built with Streamlit for rapid dashboard development
- Interactive visualizations using Plotly and Matplotlib
- Responsive UI with custom styling


## License

[MIT License](LICENSE)