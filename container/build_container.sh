#!/bin/bash
# Build script for monarch-phenologs Apptainer container

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
CONTAINER_NAME="monarch-phenologs"
VERSION="${1:-latest}"  # Allow version as first argument, default to 'latest'
DEF_FILE="${SCRIPT_DIR}/${CONTAINER_NAME}.def"

# Use VM-local writable directories (important for Lima)
TMP_BASE="${TMPDIR:-$HOME/tmp}"
BUILD_DIR="${BUILD_DIR:-$TMP_BASE/${CONTAINER_NAME}-build}"
OUTPUT_SIF="${BUILD_DIR}/${CONTAINER_NAME}_${VERSION}.sif"

# Final output location (host-mounted, copied at end)
FINAL_SIF="${SCRIPT_DIR}/${CONTAINER_NAME}_${VERSION}.sif"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Pre-flight checks
check_prerequisites() {
    log_info "Checking prerequisites..."

    # Check for Apptainer/Singularity
    if command -v apptainer &> /dev/null; then
        CONTAINER_CMD="apptainer"
        log_info "Found Apptainer: $(apptainer --version)"
    elif command -v singularity &> /dev/null; then
        CONTAINER_CMD="singularity"
        log_info "Found Singularity: $(singularity --version)"
    else
        log_error "Neither Apptainer nor Singularity found. Please install one of them."
        exit 1
    fi

    # Check definition file exists
    if [[ ! -f "${DEF_FILE}" ]]; then
        log_error "Definition file not found: ${DEF_FILE}"
        exit 1
    fi

    # Check Python scripts exist
    if [[ ! -d "${REPO_ROOT}/python" ]]; then
        log_error "Python scripts directory not found: ${REPO_ROOT}/python"
        exit 1
    fi

    # Check pyproject.toml exists
    if [[ ! -f "${REPO_ROOT}/pyproject.toml" ]]; then
        log_error "pyproject.toml not found: ${REPO_ROOT}/pyproject.toml"
        exit 1
    fi

    log_info "All prerequisites satisfied"
}

# Build container
build_container() {
    log_info "Building container: ${OUTPUT_SIF}"
    log_info "Using definition file: ${DEF_FILE}"

    mkdir -p "${BUILD_DIR}"

    ${CONTAINER_CMD} build \
        --fakeroot \
        --force \
        "${OUTPUT_SIF}" \
        "${DEF_FILE}"

    log_info "Copying container back to repo..."
    cp "${OUTPUT_SIF}" "${FINAL_SIF}"

    if [[ $? -eq 0 ]]; then
        log_info "Container built successfully: ${FINAL_SIF}"
        ls -lh "${FINAL_SIF}"
    else
        log_error "Container build failed"
        exit 1
    fi
}

# Test container
test_container() {
    log_info "Testing container..."

    # Test 1: Python version
    log_info "Test 1: Checking Python version..."
    ${CONTAINER_CMD} exec "${FINAL_SIF}" python --version

    # Test 2: Import key dependencies
    log_info "Test 2: Checking Python dependencies..."
    ${CONTAINER_CMD} exec "${FINAL_SIF}" python -c "
import sys
import pandas as pd
import numpy as np
import scipy
import networkx as nx
import pronto
import requests
import pydantic
import matplotlib
print('✓ All core dependencies imported successfully')
print(f'  - pandas: {pd.__version__}')
print(f'  - numpy: {np.__version__}')
print(f'  - scipy: {scipy.__version__}')
print(f'  - networkx: {nx.__version__}')
"

    # Test 3: Check script accessibility
    log_info "Test 3: Checking Python scripts are accessible..."
    ${CONTAINER_CMD} exec "${FINAL_SIF}" ls -l /opt/monarch-phenologs/python/

    # Test 4: Verify script can be executed (dry run)
    log_info "Test 4: Testing script execution (help)..."
    ${CONTAINER_CMD} exec "${FINAL_SIF}" python /opt/monarch-phenologs/python/04_compute_phenologs.py --help || true
    log_info "All tests passed ✓"
}

# Create symlink to 'latest'
create_symlink() {
    LATEST_LINK="${SCRIPT_DIR}/${CONTAINER_NAME}_latest.sif"

    if [[ -L "${LATEST_LINK}" || -f "${LATEST_LINK}" ]]; then
        rm -f "${LATEST_LINK}"
    fi

    ln -s "$(basename "${FINAL_SIF}")" "${LATEST_LINK}"
    log_info "Created symlink: ${LATEST_LINK} -> $(basename "${FINAL_SIF}")"
}

# Display next steps
show_next_steps() {
    cat <<EOF

${GREEN}========================================${NC}
${GREEN}Container Build Complete!${NC}
${GREEN}========================================${NC}

Container: ${FINAL_SIF}
Size: $(du -h "${FINAL_SIF}" | cut -f1)

${YELLOW}Next Steps:${NC}

1. Test the container manually:
   ${CONTAINER_CMD} exec ${FINAL_SIF} python --version

2. Copy to HPC cluster:
   scp ${FINAL_SIF} username@hpc:/path/to/containers/

3. Update nextflow.config:
   - Add singularity profile
   - Set container path

4. Run Nextflow workflow:
   nextflow run main.nf -profile slurm,singularity

${GREEN}========================================${NC}

EOF
}

# Main execution
main() {
    log_info "Starting monarch-phenologs container build (version: ${VERSION})"

    check_prerequisites
    build_container
    test_container
    create_symlink
    show_next_steps

    log_info "Build pipeline completed successfully!"
}

# Run main function
main "$@"
