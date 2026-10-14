FROM python:3.10-slim

# The locked environment, installed by uv into the image's python. Nothing is resolved at build
# time, so the image is exactly what uv.lock says
COPY --from=ghcr.io/astral-sh/uv:0.6.13 /uv /bin/uv
ENV UV_PYTHON_DOWNLOADS=never UV_LINK_MODE=copy UV_COMPILE_BYTECODE=1

WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen --all-groups

# The code. Scripts find their imports next to themselves, so they run from any working directory:
#   python /app/python/02_kg_data_extraction.py -p /data
COPY . .
ENV PATH="/app/.venv/bin:$PATH"

CMD ["/bin/bash"]
