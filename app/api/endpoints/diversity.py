from fastapi import APIRouter, HTTPException

from app.models.schemas import (
    DiversityAnalysisRequest,
    DiversityAnalysisResponse,
    MicrobiomeData,
)
from app.services.diversity_analysis import run_diversity_analysis

router = APIRouter()


@router.post("/analyze", response_model=DiversityAnalysisResponse)
async def analyze_diversity(request: DiversityAnalysisRequest, data: MicrobiomeData):
    """
    Analyze the diversity of the microbiome data

    Calculates alpha diversity metrics for each sample

    Args:
        request: parameters for diversity analysis
        data: microbiome data in standardized format

    Returns:
        Diversity analysis results
    """

    if data is None:
        raise HTTPException(
            status_code=400, detail="No data provided. Please provide data."
        )

    # Validate request parameters
    supported_metrics = ["shannon", "simpson", "pielou", "chao1"]
    for metric in request.metrics:
        if metric not in supported_metrics:
            raise HTTPException(
                status_code=400,
                detail=f"Unsupported diversity metric: {metric}.",
            )
    # Run diversity analysis
    try:
        alpha_diversity_metrics, beta_diversity_metrics, group_comparison_metrics = (
            run_diversity_analysis(data, metrics=request.metrics)
        )
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Error during diversity analysis: {str(e)}"
        )

    # Returnresponse object with all metrics
    response = DiversityAnalysisResponse(
        alpha_diversity_metrics=alpha_diversity_metrics,
        beta_diversity_metrics=beta_diversity_metrics,
        group_comparison_metrics=group_comparison_metrics,
        sample_ids=data.sample_ids,
    )

    return response
