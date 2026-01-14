#!/bin/bash

# Build script for monarch-phenologs Docker image

# Set image name and tag
IMAGE_NAME="monarch-phenologs"
TAG="latest"
USERNAME="dockerman033"

echo "Building Docker image: ${IMAGE_NAME}:${TAG}"

# Build the Docker image
##docker buildx build --platform linux/amd64,linux/arm64 -t ${USERNAME}/${IMAGE_NAME}:${TAG} .
##docker buildx build --platform linux/amd64,linux/arm64 -t ${IMAGE_NAME}:${TAG} .
docker buildx build --platform linux/amd64,linux/arm64 -t ${USERNAME}/${IMAGE_NAME}:${TAG} --push .

if [ $? -eq 0 ]; then
    echo "Docker image built successfully!"
    echo ""
    echo "To run the container:"
    echo "  docker run -it --rm ${IMAGE_NAME}:${TAG}"
    echo ""
    echo "To run with a data volume mounted:"
    echo "  docker run -it --rm -v \$(pwd)/data:/data ${IMAGE_NAME}:${TAG}"
    echo ""
    echo "To run a specific Python script:"
    echo "  docker run -it --rm -v \$(pwd)/data:/data ${IMAGE_NAME}:${TAG} python python/01_download_and_setup.py"
else
    echo "Docker build failed!"
    exit 1
fi
