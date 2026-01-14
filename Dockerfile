# Use Python 3.10 slim image as base
FROM python:3.10-slim

# Set working directory
WORKDIR /app

# # Install system dependencies that may be needed for scientific Python packages
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    && rm -rf /var/lib/apt/lists/*

# Assumes poetry is used for dependency management
# Copy only dependency files first (better caching)
COPY pyproject.toml poetry.lock* ./

# Install poetry and generate requirements.txt
RUN pip install --no-cache-dir poetry poetry-plugin-export && \
    poetry export --format requirements.txt --without-hashes -o requirements.txt

    
# Install Python dependencies and clean up
RUN pip install --no-cache-dir -r requirements.txt && \
    pip uninstall -y poetry poetry-plugin-export && \
    rm -rf ~/.cache/pip

# Copy the entire project (so whole code base is available)
COPY . .

# Set Python path to include the app directory
ENV PYTHONPATH=/app:$PYTHONPATH

# Create a directory for data/output (can be mounted as volume)
RUN mkdir -p /data

# Set default command to bash for interactive use
CMD ["/bin/bash"]
