from fastapi import FastAPI, HTTPException

from app.api.endpoints import upload

app = FastAPI(
    title="Microbiome Analysis API",
    description="API for calculating various diversity metrics in microbiome data",
    version="0.1.0"
)

app.include_router(upload.router, prefix="/upload")