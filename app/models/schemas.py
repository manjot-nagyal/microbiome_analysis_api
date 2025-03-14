from typing import Dict, List, Optional, Union

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

    metrics: List[str] = Field(
        ["shannon", "simpson", "pielou", "chao1"],
        description="List of diversity metrics to calculate",
    )

    group_by: Optional[str] = Field(
        None, description="Optional metadata column to group samples by for comparison"
    )


class DiversityAnalysisResponse(BaseModel):
    """Response schema for diversity analysis."""

    alpha_diversity_metrics: Dict[str, Dict[str, float]] = Field(
        ..., description="Alpha diversity metrics for each sample id"
    )

    beta_diversity_metrics: Optional[Dict[str, List[List[float]]]] = Field(
        ..., description="Beta diversity metrics for each sample id"
    )

    # group comparison metrics
    group_comparison_metrics: Optional[Dict[str, Dict[str, Dict[str, float]]]] = Field(
        ..., description="Group comparison metrics for each sample id"
    )

    sample_ids: List[str] = Field(..., description="List of sample ids")


class CorrelationAnalysisRequest(BaseModel):
    """Request schema for correlation analysis."""

    correlation_method: str = Field(
        "spearman",
        description="Correlation method to use ('pearson', 'spearman', or 'kendall')",
    )
    metadata_columns: Optional[List[str]] = Field(
        None, description="Optional list of metadata columns to correlate with taxa"
    )
    min_abundance: float = Field(
        0.01,
        description="Minimum relative abundance threshold for "
        "a taxon to be included in the analysis",
    )
    min_prevalence: float = Field(
        0.1,
        description="Minimum fraction of samples in which a "
        "taxon must be present to be included",
    )


class CorrelationAnalysisResponse(BaseModel):
    """Response schema for correlation analysis results."""

    taxon_correlations: Dict[str, Dict[str, float]] = Field(
        ..., description="Correlations between taxa, keyed by pairs of taxonomy IDs"
    )
    metadata_correlations: Optional[Dict[str, Dict[str, float]]] = Field(
        None,
        description="Correlations between taxa and metadata, "
        "keyed by taxonomy ID and metadata column",
    )
    p_values: Dict[str, Dict[str, float]] = Field(
        ...,
        description="P-values for correlations, with the same "
        "structure as the correlation matrices",
    )
    filtered_taxa: List[str] = Field(
        ...,
        description="List of taxonomy IDs that passed the "
        "abundance and prevalence filters",
    )


class AnalysisResult(BaseModel):
    """Generic schema for analysis results."""

    analysis_type: str
    result: Union[DiversityAnalysisResponse, CorrelationAnalysisResponse]
