# Fragment Shift Estimation Pipeline - Complete Setup and Execution Guide

This comprehensive guide provides step-by-step instructions for setting up and running the cross-correlation-based fragment shift estimation pipeline in WSL.

## Overview

This pipeline provides:
- **Cross-correlation analysis** using FFT acceleration for performance
- **Robust fragment shift estimation** with confidence scoring
- **Shift-aware peak calling** for improved accuracy
- **Comprehensive validation** with synthetic test datasets
- **Multi-sample comparison** and quality assessment

## Prerequisites

- **WSL Environment**: Ubuntu 20.04+ or similar Linux distribution
- **System Resources**: 8GB RAM, 15GB free disk space
- **Network**: Internet connection for dependencies
- **Time**: 60-90 minutes for complete setup and testing

## Step 1: System Environment Setup

### Install Core Dependencies

```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install essential build tools
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
    python3-venv \
    python3-numpy \
    python3-scipy \
    openjdk-11-jdk \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfftw3-dev \
    bc \
    jq

# Install additional Python packages
pip3 install --user numpy scipy matplotlib seaborn pandas
```

### Install Rust with Optimization Features

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source ~/.cargo/env

# Verify installation
rustc --version
cargo --version

# Install additional Rust components for optimization
rustup component add llvm-tools-preview
```

### Install Bioinformatics Tools

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda3
export PATH="$HOME/miniconda3/bin:$PATH"
echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Initialize conda
conda init bash
source ~/.bashrc

# Create bioinformatics environment
conda create -n biotools -y python=3.9
conda activate biotools

# Install bioinformatics tools
conda install -y -c bioconda \
    fastqc=0.11.9 \
    fastp=0.23.2 \
    bowtie2=2.5.1 \
    samtools=1.17 \
    bedtools=2.30.0

# Install Nextflow
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
sudo chmod +x /usr/local/bin/nextflow

# Verify installations
nextflow -version
bowtie2 --version
samtools --version
```

## Step 2: Project Setup

### Create Project Structure

```bash
# Create main project directory
mkdir -p ~/fragment_shift_pipeline
cd ~/fragment_shift_pipeline

# Initialize Rust project
cargo init --name fragment-shift-estimator

# Create additional directories
mkdir -p {src,docs,examples,benchmarks}
```

### Set Up Project Files

#### 1. Copy Cargo.toml

```bash
# Copy the Cargo.toml content from the artifact
cat > Cargo.toml << 'EOF'
[package]
name = "fragment-shift-estimator"
version = "0.4.0"
edition = "2021"
authors = ["ChIP-seq Analysis Team"]
description = "Cross-correlation based fragment shift estimation for ChIP-seq data"

[[bin]]
name = "shift_estimator"
path = "src/main.rs"

[[bin]]
name = "peak_caller"
path = "src/peak_caller.rs"

[[bin]]
name = "shift_reads"
path = "src/shift_reads.rs"

[[bin]]
name = "coverage_analyzer"
path = "src/coverage_analyzer.rs"

[dependencies]
rust-htslib = "0.46"
rayon = "1.8"
serde = { version = "1.0", features = ["derive"] }
serde_json = "1.0"
clap = { version = "4.0", features = ["derive"] }
anyhow = "1.0"
log = "0.4"
env_logger = "0.10"
bio = "1.6"
itertools = "0.12"
statrs = "0.16"
ndarray = "0.15"
dashmap = "5.5"
crossbeam = "0.8"
indicatif = "0.17"
chrono = { version = "0.4", features = ["serde"] }
flate2 = "1.0"
csv = "1.3"
regex = "1.10"
num-complex = "0.4"
rustfft = "6.1"

[dev-dependencies]
tempfile = "3.0"
criterion = "0.5"
approx = "0.5"

[[bench]]
name = "shift_benchmarks"
harness = false

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
panic = "abort"
debug = false

[profile.dev]
opt-level = 1
debug = true

[features]
default = ["fft-acceleration"]
fft-acceleration = ["rustfft"]
simd-optimization = []
EOF
```

#### 2. Copy Source Files

```bash
# Copy main.rs to src/main.rs (from artifact)
# Copy peak_caller.rs to src/peak_caller.rs (from artifact)
# Copy main.nf (from artifact)
# Copy generate_shift_test_data.py (from artifact)

# Make scripts executable
chmod +x generate_shift_test_data.py
```

## Step 3: Build and Test Rust Components

### Compile with Optimizations

```bash
# Navigate to project directory
cd ~/fragment_shift_pipeline

# Clean any previous builds
cargo clean

# Build with release optimizations (includes FFT)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# This compilation will take 10-20 minutes on first run
# Progress will be shown for each dependency

# Verify binaries are created
ls -la target/release/
file target/release/shift_estimator
```

### Test Basic Functionality

```bash
# Test shift estimator help
cargo run --release --bin shift_estimator -- --help

# Test peak caller help  
cargo run --release --bin peak_caller -- --help

# Run unit tests
cargo test --release

# Run with verbose output if needed
cargo test --release -- --nocapture
```

## Step 4: Generate Comprehensive Test Dataset

### Create Fragment Shift Test Data

```bash
# Generate test dataset with multiple known fragment shifts
python3 generate_shift_test_data.py \
    --output-dir shift_test \
    --num-reads 150000 \
    --read-length 75 \
    --fragment-shifts 100 150 200 250 \
    --fragment-size-std 25 \
    --binding-sites-per-mb 20.0 \
    --seed 42

# Verify generated files
ls -la shift_test/
ls -la shift_test/data/
ls -la shift_test/reference/
```

### Inspect Generated Test Data

```bash
cd shift_test

# Check reference genome
head reference/genome.fa
samtools faidx reference/genome.fa
head reference/genome.fa.fai

# Check binding sites
head -10 reference/binding_sites.bed
wc -l reference/binding_sites.bed

# Check sample descriptions
cat data/*_description.txt

# Verify FASTQ files
zcat data/sample_shift_150bp.fastq.gz | head -8

# Check file sizes
du -sh data/*.fastq.gz
```

## Step 5: Test Individual Pipeline Components

### Test Reference Indexing and Alignment

```bash
# Activate conda environment
conda activate biotools

# Build reference indices
echo "Building bowtie2 index..."
bowtie2-build reference/genome.fa reference/genome_index

echo "Building samtools index..."
samtools faidx reference/genome.fa

# Test alignment with one sample
echo "Testing alignment pipeline..."
sample="sample_shift_150bp"

# Quality control and trimming
fastp \
    -i data/${sample}.fastq.gz \
    -o ${sample}_trimmed.fastq.gz \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --thread 4 \
    --json ${sample}_fastp.json \
    --html ${sample}_fastp.html

# Alignment
bowtie2 \
    -x reference/genome_index \
    -U ${sample}_trimmed.fastq.gz \
    --threads 4 \
    --very-sensitive \
    | samtools view -bS -F 4 -q 10 - \
    | samtools sort -@ 4 -o ${sample}_aligned.bam -

# Index BAM
samtools index ${sample}_aligned.bam

# Check alignment statistics
samtools flagstat ${sample}_aligned.bam
samtools stats ${sample}_aligned.bam | head -20
```

### Test Fragment Shift Estimation

```bash
# Test shift estimator on aligned sample
echo "Testing fragment shift estimation..."
cargo run --release --bin shift_estimator -- \
    --input ${sample}_aligned.bam \
    --output ${sample}_shift_estimate.json \
    --correlation-output ${sample}_correlation_profile.json \
    --max-shift 400 \
    --bin-size 1 \
    --use-fft \
    --smoothing-window 5 \
    --threads 4 \
    --verbose

# Check results
cat ${sample}_shift_estimate.json | jq .
echo "Expected shift: 150 bp"
echo "Estimated shift: $(jq -r '.estimated_shift' ${sample}_shift_estimate.json) bp"
echo "Confidence: $(jq -r '.confidence_score' ${sample}_shift_estimate.json)"
```

### Test Shift-Aware Peak Calling

```bash
# Test peak calling with shift correction
echo "Testing shift-aware peak calling..."
cargo run --release --bin peak_caller -- \
    ${sample}_aligned.bam \
    ${sample}_shift_estimate.json \
    --output ${sample}_peaks.bed \
    --window-size 200 \
    --min-coverage 3.0 \
    --pvalue-threshold 0.05 \
    --threads 4 \
    --verbose

# Check peak results
head ${sample}_peaks.bed
wc -l ${sample}_peaks.bed
```

## Step 6: Run Complete Pipeline

### Execute Full Nextflow Pipeline

```bash
# Navigate to test directory
cd shift_test

# Run complete pipeline
echo "Running complete fragment shift estimation pipeline..."
nextflow run ../main.nf \
    -c nextflow.config \
    --input_dir data \
    --output_dir results \
    --reference_genome reference/genome.fa \
    --max_shift 500 \
    --use_fft true \
    --bin_size 1 \
    --threads 4 \
    -with-report results/execution_report.html \
    -with-timeline results/timeline.html