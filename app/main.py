from fastapi import FastAPI

from app.api.endpoints import diversity, upload

app = FastAPI(
    title="Microbiome Analysis API",
    description="API for microbiome data analysis",
    version="0.1.0",
)

app.include_router(upload.router, prefix="/upload")
app.include_router(diversity.router, prefix="/diversity")
