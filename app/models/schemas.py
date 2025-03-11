from typing import List, Optional, Dict, Any

from pydantic import BaseModel, Field


# Create a MicrobiomeData class
class MicrobiomeData(BaseModel):
    """
    A Schema for the MicrobiomeData object
    """

    sample_ids: List[str] = Field(
        ..., description="List of sample identifiers")
    
    feature_ids: List[str] = Field(
        ..., description="List of feature identifiers (i.e. OTUs, ASVs, etc.)")
    
    counts_matrix: List[List[float]] = Field(
        ..., description="Counts matrix with samples as rows and features as columns")
    
    metadata: Optional[Dict[str, Any]] = Field(
        None, description="Optional metadata for samples, where keys are sample_ids and values are dictionaries of metadata")
