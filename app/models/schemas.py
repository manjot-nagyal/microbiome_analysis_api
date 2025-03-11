from typing import List, Optional, Dict, Union

from pydantic import BaseModel, Field


# Create a MicrobiomeData class
class MicrobiomeData(BaseModel):
    """
    A Schema for the MicrobiomeData object
    """

    sample_ids: List[str] = Field(..., description="List of sample ids")

    feature_ids: List[str] = Field(
        ..., description="List of feature ids (i.e. OTUs, ASVs, etc.)"
    )

    counts_matrix: List[List[float]] = Field(
        ...,
        description="Counts matrix (samples by features)",
    )

    metadata: Optional[Dict[str, Dict[str, Union[str, float, int, bool]]]] = Field(
        None,
        description="Optional metadata for samples, "
        "where keys are sample_ids and values are dictionaries "
        "of metadata",
    )


class DiversityAnalysisRequest(BaseModel):
    """Request schema for diversity analysis."""

    diversity_metrics: List[str] = Field(
        ["shannon", "simpson", "pielou", "chao1"],
        description="List of diversity metrics to calculate",
    )
