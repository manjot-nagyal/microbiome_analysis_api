FROM python:3.12-slim AS base

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=off \
    PIP_DISABLE_PIP_VERSION_CHECK=on

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    gcc \
    libc6-dev \
    libhdf5-dev \
    python3-dev \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install uv


# Test
FROM base AS test

COPY . .
RUN uv pip install --system -e ".[dev]"

ARG TEST_CMD="pytest -v"
RUN ${TEST_CMD}


# Backend
FROM base AS backend

COPY . .
RUN uv pip install --system -e "."

RUN apt-get update && apt-get install -y curl && apt-get clean && rm -rf /var/lib/apt/lists/*

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]


# Frontend
FROM base AS frontend

COPY . .
RUN uv pip install --system -e "."

EXPOSE 8501

CMD uv run streamlit run frontend.py -- --server.port=8501 --server.address=0.0.0.0
