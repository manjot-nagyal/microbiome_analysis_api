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

COPY . .
RUN uv pip install --system -e "."

# Test
FROM base AS test

RUN uv pip install --system -e ".[all]"

ARG TEST_CMD="pytest -v"
RUN ${TEST_CMD}


# Backend
FROM base AS backend

RUN apt-get update && apt-get install -y curl && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN uv pip install --system -e ".[backend]"

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]


# Frontend
FROM base AS frontend

RUN uv pip install --system -e ".[frontend]"

EXPOSE 8501

CMD ["streamlit", "run", "frontend.py", "--", "--server.port=8501", "--server.address=0.0.0.0"]
