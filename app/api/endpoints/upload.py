import json
from typing import Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from app.models.schemas import MicrobiomeData
from app.utils.parsers import normalize_counts, parse_csv_data, parse_json_data

router = APIRouter()


current_data = None


@router.post("/", response_model=MicrobiomeData)
async def upload_microbiome_data(
    file: Optional[UploadFile] = File(None),
    json_data: Optional[str] = Form(None),
    normalize: bool = Form(True),
):
    """
    Upload microbiome data as CSV format or JSON

    Returns parsed microbiome data in a standardized format

    ## File format requirements:

    ### CSV structure
    The CSV should have the following structure:
    1st column: sample_id
    1st row: feature_id (e.g. OTU, ASV, species names...)
    Remaining cells: feature_counts
    Optional: metadata columsn with the 'metadata_' prefix

    ### Requirements
    - Each sample must have abundance values for all taxa
    - The dimensions of the abundance matrix must match the number of
    samples and taxa
    - Metadata is optional but if provided must be keyed by sample IDs
    """
    global current_data

    if file is None and json_data is None:
        raise HTTPException(
            status_code=400,
            detail="No file uploaded or JSON data provided",
        )

    try:
        if file is not None:
            try:
                content = await file.read()
                content_str = content.decode("utf-8")
            except UnicodeDecodeError:
                raise HTTPException(
                    status_code=400,
                    detail="File encoding error. "
                    "File is not a valid CSV. "
                    "Please ensure the file is in UTF-8 format.",
                )

            # For CSV files
            if file.filename and file.filename.endswith(".csv"):
                try:
                    print(content_str)
                    data = parse_csv_data(content_str)

                except Exception as e:

                    raise HTTPException(
                        status_code=400,
                        detail=f"Error parsing CSV file: {str(e)}. "
                        "Please check the file format and structure "
                        "against the documentation.",
                    )
            # For JSON files
            elif file.filename and file.filename.endswith(".json"):
                try:
                    data = parse_json_data(content_str)
                except json.JSONDecodeError:
                    raise HTTPException(
                        status_code=400,
                        detail="Invalid JSON format. Please ensure "
                        "the file contains valid JSON data.",
                    )
                except Exception as e:
                    raise HTTPException(
                        status_code=400,
                        detail=f"Error parsing JSON data: {str(e)}. Please check "
                        "file format/structure against documentation.",
                    )

            else:
                raise HTTPException(
                    status_code=400,
                    detail="Unsupported file type. "
                    "Please upload a correct file format.",
                )
        # Directly from JSON data
        else:
            try:
                if json_data is None:
                    raise HTTPException(
                        status_code=400,
                        detail="No JSON data provided. Please provide valid JSON data.",
                    )
                data = parse_json_data(json_data)
            except json.JSONDecodeError:
                raise HTTPException(
                    status_code=400,
                    detail="Invalid JSON data. Please ensure you are "
                    "submitting valid JSON in the required format.",
                )
            except Exception as e:
                raise HTTPException(
                    status_code=400,
                    detail=f"Error parsing JSON data: {str(e)}. Please check "
                    "your data against the required format in the documentation.",
                )

        if normalize:
            data = normalize_counts(data)

        current_data = data
        return data
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error parsing data: {str(e)}")


@router.post("/current")
async def reset_current_data():
    """
    Reset the currently loaded microbiome data.

    This endpoint is useful for clearing the current data between analyses.
    """
    global current_data
    current_data = None
    return {"message": "Current microbiome data has been reset."}


@router.get("/current", response_model=Optional[MicrobiomeData])
async def get_current_data():
    """
    Get the currently loaded microbiome data.

    Returns None if no data has been uploaded yet.
    """
    global current_data
    return current_data
