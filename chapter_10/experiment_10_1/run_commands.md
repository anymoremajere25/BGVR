# Epigenomic Peak Calling Pipeline - Execution Guide

This guide provides step-by-step commands to set up and run the epigenomic peak calling pipeline in WSL.

## Prerequisites

- WSL (Windows Subsystem for Linux) with Ubuntu 20.04+ or similar Linux distribution
- At least 8GB RAM and 10GB free disk space
- Internet connection for downloading dependencies

## Step 1: Initial Setup

```bash
# Update your system
sudo apt update && sudo apt upgrade -y

# Create project directory
mkdir -p ~/epigenomic_pipeline
cd ~/epigenomic_pipeline

# Make setup script executable and run it
chmod +x setup_environment.sh
./setup_environment.sh
```

## Step 2: Set Up Project Files

After running the setup script, copy the provided files:

```bash
# Navigate to project directory
cd ~/epigenomic_pipeline

# Copy Cargo.toml (create this file with the content from artifact)
cat > Cargo.toml << 'EOF'
[package]
name = "epigenomic-peak-caller"
version = "0.1.0"
edition = "2021"

[[bin]]
name = "peak_caller"
path = "src/main.rs"

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

[dev-dependencies]
tempfile = "3.0"
EOF

# Create src directory and copy main.rs
mkdir -p src
# Copy the main.rs content from the artifact to src/main.rs

# Copy main.nf (Nextflow pipeline file)
# Copy the main.nf content from the artifact

# Copy generate_test_data.py
# Copy the generate_test_data.py content from the artifact
chmod +x generate_test_data.py
```

## Step 3: Build Rust Dependencies

```bash
# Restart terminal or source environment
source ~/.bashrc

# Activate conda environment
conda activate biotools

# Build Rust project (this will download and compile dependencies)
cargo build --release

# Test the Rust binary
cargo run --bin peak_caller -- --help
```

## Step 4: Generate Test Dataset

```bash
# Generate synthetic test data
python3 generate_test_data.py \
    --output-dir epigenomic_test \
    --num-reads 100000 \
    --num-samples 3 \
    --read-length 75 \
    --seed 42

# Check generated files
ls -la epigenomic_test/
ls -la epigenomic_test/data/
ls -la epigenomic_test/reference/
```

## Step 5: Test Individual Components

### Test Rust Peak Caller Directly

```bash
# Navigate to test data directory
cd epigenomic_test

# First, we need to create a BAM file from FASTQ for testing
# Build bowtie2 index
bowtie2-build reference/genome.fa reference/genome_index

# Align one sample to create test BAM
fastp -i data/sample_01.fastq.gz -o sample_01_trimmed.fastq.gz
bowtie2 -x reference/genome_index -U sample_01_trimmed.fastq.gz | \
    samtools view -bS - | \
    samtools sort -o sample_01_aligned.bam -
samtools index sample_01_aligned.bam

# Test Rust peak caller
cargo run --release --bin peak_caller -- \
    --input sample_01_aligned.bam \
    --output test_peaks.json \
    --window-size 200 \
    --min-coverage 3.0 \
    --pvalue-threshold 0.05 \
    --threads 4

# Check results
head test_peaks.json
```

### Test Nextflow Pipeline

```bash
# Run the complete Nextflow pipeline
cd ~/epigenomic_pipeline/epigenomic_test

# Run pipeline with test configuration
nextflow run ../main.nf \
    --input_dir data \
    --output_dir results \
    --reference_genome reference/genome.fa \
    --window_size 200 \
    --min_coverage 3.0 \
    --pvalue_threshold 0.05 \
    --threads 4

# Check pipeline results
ls -la results/
```

## Step 6: Examine Results

```bash
# View pipeline outputs
cd results

# Check QC results
ls qc/fastqc/
ls qc/trimmed/

# Check alignments
ls alignments/
samtools flagstat alignments/sample_01_aligned.bam

# Check peak calling results
ls peaks/
head peaks/sample_01_peaks.json
head peaks/sample_01_peaks.bed

# Check merged results
ls merged/
head merged/merged_peaks.bed

# View HTML report
# If you have a web browser available:
# firefox pipeline_report.html
# Or copy the file to Windows and open it there
```

## Step 7: Advanced Usage Examples

### Run with Custom Parameters

```bash
# Run with more stringent parameters
nextflow run ../main.nf \
    --input_dir data \
    --output_dir results_stringent \
    --reference_genome reference/genome.fa \
    --window_size 150 \
    --min_coverage 10.0 \
    --pvalue_threshold 0.01 \
    --threads 8
```

### Run Rust Peak Caller with Different Settings

```bash
# High sensitivity (more peaks)
cargo run --release --bin peak_caller -- \
    --input sample_01_aligned.bam \
    --output peaks_sensitive.json \
    --window-size 100 \
    --min-coverage 2.0 \
    --pvalue-threshold 0.1 \
    --threads 4

# High specificity (fewer, more confident peaks)
cargo run --release --bin peak_caller -- \
    --input sample_01_aligned.bam \
    --output peaks_specific.json \
    --window-size 300 \
    --min-coverage 10.0 \
    --pvalue-threshold 0.001 \
    --threads 4
```

## Step 8: Performance Monitoring

```bash
# Monitor resource usage during pipeline execution
# In another terminal:
htop

# Check disk usage
df -h
du -sh epigenomic_test/

# Monitor log files
tail -f .nextflow.log
```

## Troubleshooting

### Common Issues and Solutions

1. **Rust compilation errors:**
   ```bash
   # Update Rust toolchain
   rustup update
   
   # Clean and rebuild
   cargo clean
   cargo build --release
   ```

2. **Conda environment issues:**
   ```bash
   # Recreate conda environment
   conda deactivate
   conda env remove -n biotools
   conda create -n biotools python=3.9
   conda activate biotools
   conda install -y -c bioconda fastqc fastp bowtie2 samtools bedtools
   ```

3. **Memory issues:**
   ```bash
   # Reduce number of reads in test data
   python3 generate_test_data.py --num-reads 50000 --num-samples 2
   
   # Reduce number of threads
   cargo run --bin peak_caller -- --threads 2
   ```

4. **Permission errors:**
   ```bash
   # Fix file permissions
   chmod +x *.sh
   chmod +x *.py
   sudo chown -R $USER:$USER ~/epigenomic_