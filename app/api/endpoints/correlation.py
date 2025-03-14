from fastapi import APIRouter, Depends, HTTPException

from app.api.endpoints.upload import get_current_data
from app.models.schemas import CorrelationAnalysisRequest, CorrelationAnalysisResponse
from app.services.correlation_analysis import run_correlation_analysis

router = APIRouter(prefix="/correlation")


@router.post("/analyze", response_model=CorrelationAnalysisResponse)
async def analyze_correlation(
    request: CorrelationAnalysisRequest, data=Depends(get_current_data)
):
    """
    Analyze correlations in microbiome data.

    Calculates correlations between different taxa and between taxa and metadata
    variables if specified. Taxa with low abundance or prevalence are filtered
    out to focus on relevant relationships.

    Args:
        request: Parameters for correlation analysis
        data: Microbiome data from the upload endpoint

    Returns:
        Correlation analysis results
    """
    if data is None:
        raise HTTPException(
            status_code=400,
            detail="No microbiome data available. Please upload data first.",
        )

    supported_methods = ["pearson", "spearman", "kendall"]
    if request.correlation_method not in supported_methods:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported correlation method: {request.correlation_method}. "
            "Supported methods: {', '.join(supported_methods)}",
        )

    if request.metadata_columns:
        if not data.metadata:
            raise HTTPException(
                status_code=400,
                detail="No metadata available in the uploaded data, but metadata "
                "columns were requested for correlation.",
            )

        available_columns = set()
        for sample_id, metadata in data.metadata.items():
            available_columns.update(metadata.keys())

        unknown_columns = set(request.metadata_columns) - available_columns
        if unknown_columns:
            raise HTTPException(
                status_code=400,
                detail=f"Unknown metadata columns: {', '.join(unknown_columns)}. "
                "Available columns: {', '.join(available_columns)}",
            )

    try:
        taxon_correlations, metadata_correlations, p_values, filtered_taxa = (
            run_correlation_analysis(
                data=data,
                correlation_method=request.correlation_method,
                metadata_columns=request.metadata_columns,
                min_abundance=request.min_abundance,
                min_prevalence=request.min_prevalence,
            )
        )

        if not filtered_taxa:
            raise HTTPException(
                status_code=400,
                detail="No taxa passed the abundance and prevalence filters. "
                "Try lowering the thresholds.",
            )

        response = CorrelationAnalysisResponse(
            taxon_correlations=taxon_correlations,
            metadata_correlations=metadata_correlations,
            p_values=p_values,
            filtered_taxa=filtered_taxa,
        )

        return response

    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Error in correlation analysis: {str(e)}"
        )
