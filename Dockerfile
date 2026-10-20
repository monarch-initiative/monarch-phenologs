FROM python:3.10-slim

# ps, which Nextflow's task wrapper needs
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*

# The locked environment, installed by uv into the image's python. Nothing is resolved at build
# time, so the image is exactly what uv.lock says. uv's download cache is not kept (it would add ~500 MB)
COPY --from=ghcr.io/astral-sh/uv:0.6.13 /uv /bin/uv
ENV UV_PYTHON_DOWNLOADS=never UV_LINK_MODE=copy UV_COMPILE_BYTECODE=1

WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --all-groups --no-cache

# The code. Scripts find their imports next to themselves, so they run from any working directory:
#   python /app/python/02_kg_data_extraction.py -p /data
COPY . .
ENV PATH="/app/.venv/bin:$PATH"

# The image's date tag and the commit it is built from, which 07_build_release.py records in each release. The Build
# Image workflow sets them; a local build can too:
#   docker build --build-arg PHENOLOGS_IMAGE=local --build-arg PHENOLOGS_COMMIT=$(git rev-parse --short HEAD) .
ARG PHENOLOGS_IMAGE=
ARG PHENOLOGS_COMMIT=
ENV PHENOLOGS_IMAGE=$PHENOLOGS_IMAGE PHENOLOGS_COMMIT=$PHENOLOGS_COMMIT

CMD ["/bin/bash"]
