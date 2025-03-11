from fastapi import APIRouter, File, UploadFile, HTTPException
from typing import Optional

import logging

logger = logging.getLogger(__name__)

from app.models.schemas import MicrobiomeData
from app.utils.parsers import parse_csv_data

router = APIRouter()


current_data = None

@router.post("/", response_model=MicrobiomeData)
async def upload_microbiome_data(file: Optional[UploadFile] = File(None)):
    """
    Upload microbiome data as CSV format

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
    - The dimensions of the abundance matrix must match the number of samples and taxa
    - Metadata is optional but if provided must be keyed by sample IDs
    """
    global current_data

    if file is None:
        raise HTTPException(
            status_code=400,
            detail="No file uploaded",
        )
    
    try:
        if file is not None:
            try:    
                content_str = await file.read()
                content_str = content_str.decode("utf-8")
            except UnicodeDecodeError:
                raise HTTPException(
                    status_code=400,
                    detail="File encoding error. File is not a valid CSV. Please ensure the file is in UTF-8 format.",
                )

            # Opens file if it has a .csv extension otherwise raises an error
            if file.filename and file.filename.endswith(".csv"):
                try:
                    data = parse_csv_data(content_str)
                except Exception as e:
                    logger.error(f"Error parsing CSV file: {str(e)}")
                    raise HTTPException(
                        status_code=400,
                        detail=f"Error parsing CSV file: {str(e)}. Please check the file format and structure against the documentation.",
                    )
            else:
                logger.error("Unsupported file type uploaded")
                raise HTTPException(
                        status_code=400,
                        detail="Unsupported file type. Please upload a correct file format.",
                    )

        current_data = data

        return data
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error parsing data: {str(e)}")
