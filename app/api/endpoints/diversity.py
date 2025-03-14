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

    print(
        f"Diversity analysis request received. "
        f"Metrics: {request.metrics}, "
        f"Group by: {request.group_by}"
    )
    print(
        f"Data received: {len(data.sample_ids)} "
        f"samples, {len(data.feature_ids)} features"
    )

    if data is None or not data.sample_ids or len(data.sample_ids) == 0:
        print("Error: No data provided")
        raise HTTPException(
            status_code=400, detail="No data provided. Please provide data."
        )

    supported_metrics = ["shannon", "simpson", "pielou", "chao1"]
    for metric in request.metrics:
        if metric not in supported_metrics:
            raise HTTPException(
                status_code=400,
                detail=f"Unsupported diversity metric: {metric}.",
            )

    try:
        print("Running diversity analysis...")
        alpha_diversity_metrics, beta_diversity_metrics, group_comparison_metrics = (
            run_diversity_analysis(data, metrics=request.metrics)
        )

        print("Diversity analysis completed successfully")
        if alpha_diversity_metrics and len(alpha_diversity_metrics) > 0:
            first_metric = list(alpha_diversity_metrics.keys())[0]
            print(
                f"Sample of {first_metric} results: "
                f"{list(alpha_diversity_metrics[first_metric].items())[:2]}"
            )

        response = DiversityAnalysisResponse(
            alpha_diversity_metrics=alpha_diversity_metrics,
            beta_diversity_metrics=beta_diversity_metrics,
            group_comparison_metrics=group_comparison_metrics,
            sample_ids=data.sample_ids,
        )

        return response

    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Error during diversity analysis: {str(e)}"
        )
