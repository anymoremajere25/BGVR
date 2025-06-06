#!/bin/bash

# Epigenomic Peak Calling Pipeline Setup Script for WSL
# This script installs all necessary dependencies and sets up the environment

set -e  # Exit on any error

echo "=========================================="
echo "Epigenomic Pipeline Environment Setup"
echo "=========================================="

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if running in WSL
if grep -q microsoft /proc/version; then
    print_status "Detected WSL environment"
else
    print_warning "This script is optimized for WSL, but will attempt to run anyway"
fi

# Update system packages
print_status "Updating system packages..."
sudo apt update && sudo apt upgrade -y

# Install essential build tools
print_status "Installing build essentials..."
sudo apt install -y \
    build-essential \
    cmake \
    pkg-config \
    git \
    wget \
    curl \
    unzip \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    openjdk-11-jdk \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev

# Install Rust
if ! command -v rustc &> /dev/null; then
    print_status "Installing Rust..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source ~/.cargo/env
    rustup update
else
    print_status "Rust already installed"
fi

# Ensure Rust is in PATH
export PATH="$HOME/.cargo/bin:$PATH"
source ~/.cargo/env

# Install bioinformatics tools
print_status "Installing bioinformatics tools..."

# Install conda/miniconda if not present
if ! command -v conda &> /dev/null; then
    print_status "Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda3
    export PATH="$HOME/miniconda3/bin:$PATH"
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    rm miniconda.sh
else
    print_status "Conda already installed"
fi

# Initialize conda
export PATH="$HOME/miniconda3/bin:$PATH"
conda init bash

# Create conda environment for bioinformatics tools
print_status "Creating bioinformatics conda environment..."
conda create -n biotools -y python=3.9
conda activate biotools

# Install bioinformatics tools via conda
print_status "Installing bioinformatics software..."
conda install -y -c bioconda \
    fastqc \
    fastp \
    bowtie2 \
    samtools \
    bedtools

# Install Nextflow
if ! command -v nextflow &> /dev/null; then
    print_status "Installing Nextflow..."
    wget -qO- https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    sudo chmod +x /usr/local/bin/nextflow
else
    print_status "Nextflow already installed"
fi

# Install Python packages
print_status "Installing Python packages..."
pip3 install --user \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    biopython \
    pysam

# Create project directory structure
PROJECT_DIR="$HOME/epigenomic_pipeline"
if [ ! -d "$PROJECT_DIR" ]; then
    print_status "Creating project directory at $PROJECT_DIR"
    mkdir -p "$PROJECT_DIR"
    cd "$PROJECT_DIR"
    
    # Initialize Rust project
    print_status "Initializing Rust project..."
    cargo init --name epigenomic-peak-caller
    
    # Create necessary subdirectories
    mkdir -p {data,reference,results,scripts}
    
else
    print_status "Project directory already exists at $PROJECT_DIR"
    cd "$PROJECT_DIR"
fi

# Add environment variables to bashrc
print_status "Setting up environment variables..."
cat >> ~/.bashrc << 'EOF'

# Epigenomic Pipeline Environment
export EPIGENOMIC_PROJECT="$HOME/epigenomic_pipeline"
export PATH="$HOME/.cargo/bin:$PATH"
export PATH="$HOME/miniconda3/bin:$PATH"

# Convenient aliases
alias activate-biotools='conda activate biotools'
alias epi-cd='cd $EPIGENOMIC_PROJECT'

EOF

# Create a simple test script
cat > test_installation.sh << 'EOF'
#!/bin/bash

echo "Testing installation..."

# Test Rust
echo "Testing Rust..."
rustc --version || echo "ERROR: Rust not found"

# Test bioinformatics tools
echo "Testing bioinformatics tools..."
source ~/.bashrc
conda activate biotools

fastqc --version || echo "ERROR: FastQC not found"
fastp --version || echo "ERROR: fastp not found"
bowtie2 --version || echo "ERROR: bowtie2 not found"
samtools --version || echo "ERROR: samtools not found"
bedtools --version || echo "ERROR: bedtools not found"

# Test Nextflow
echo "Testing Nextflow..."
nextflow -version || echo "ERROR: Nextflow not found"

# Test Python packages
echo "Testing Python packages..."
python3 -c "import numpy, pandas, matplotlib, Bio, pysam; print('Python packages OK')" || echo "ERROR: Python packages missing"

echo "Installation test complete!"
EOF

chmod +x test_installation.sh

# Create a simple run script
cat > run_pipeline.sh << 'EOF'
#!/bin/bash

# Activate conda environment
source ~/.bashrc
conda activate biotools

# Set project directory
cd $HOME/epigenomic_pipeline

echo "Running epigenomic peak calling pipeline..."
echo "Current directory: $(pwd)"

# Check if test data exists
if [ ! -d "epigenomic_test" ]; then
    echo "Generating test data..."
    python3 generate_test_data.py --output-dir epigenomic_test --num-reads 100000 --num-samples 2
fi

# Build Rust project
echo "Building Rust peak caller..."
cargo build --release

# Run Nextflow pipeline
echo "Running Nextflow pipeline..."
cd epigenomic_test
nextflow run ../main.nf -c pipeline.config --input_dir data --output_dir results --reference_genome reference/genome.fa

echo "Pipeline complete! Check results in epigenomic_test/results/"
EOF

chmod +x run_pipeline.sh

# Download example files to project directory
print_status "Setting up project files..."

# We'll need to copy the files created in the artifacts
print_status "Remember to copy the following files to $PROJECT_DIR:"
echo "  - Cargo.toml"
echo "  - src/main.rs"
echo "  - main.nf"
echo "  - generate_test_data.py"

# Final setup message
cat << 'EOF'

========================================
Setup Complete!
========================================

To get started:

1. Restart your terminal or run: source ~/.bashrc

2. Navigate to project directory:
   cd ~/epigenomic_pipeline

3. Copy the pipeline files:
   - Copy Cargo.toml to the project root
   - Copy main.rs to src/main.rs
   - Copy main.nf to the project root
   - Copy generate_test_data.py to the project root

4. Test the installation:
   ./test_installation.sh

5. Generate test data and run pipeline:
   ./run_pipeline.sh

Environment Variables:
- EPIGENOMIC_PROJECT: ~/epigenomic_pipeline
- Conda environment 'biotools' contains bioinformatics tools

Useful commands:
- activate-biotools: Activate conda environment
- epi-cd: Go to project directory

EOF

print_status "Setup script completed successfully!"