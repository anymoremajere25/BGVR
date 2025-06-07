## Experiment_10_6 Advanced Epigenomic Topics
## 🧬 Epigenomic HMM Segmentation Pipeline

A robust and comprehensive Hidden Markov Model (HMM) implementation in Rust for epigenomic chromatin state segmentation, with synthetic data generation and evaluation tools.

## 🚀 Quick Start

```bash
# 1. Setup (run once)
chmod +x setup_and_run.sh
./setup_and_run.sh setup --with-nextflow

# 2. Run complete pipeline
./setup_and_run.sh run-pipeline

# 3. View results
./setup_and_run.sh show-results
```

## 📋 Overview

This pipeline implements a sophisticated HMM for segmenting epigenomic tracks into chromatin states. It includes:

- **Rust Implementation**: High-performance HMM with EM algorithm training
- **Synthetic Data Generation**: Realistic epigenomic data with ground truth
- **Comprehensive Evaluation**: Multiple metrics and visualizations
- **Nextflow Integration**: Scalable workflow management
- **WSL Support**: Optimized for Windows Subsystem for Linux

## 🏗️ Architecture

```
📁 Project Structure
├── src/
│   └── main.rs                    # Main Rust HMM implementation
├── Cargo.toml                     # Rust dependencies
├── main.nf                        # Nextflow pipeline
├── generate_data.py               # Synthetic data generator
├── evaluate_segmentation.py       # Results evaluation
├── generate_pipeline_report.py    # Report generator
├── setup_and_run.sh              # Setup and execution script
└── README.md                      # This file
```

## 🔧 Installation

### Prerequisites

- Ubuntu/WSL2 environment
- Internet connection for downloading dependencies

### Automatic Setup

```bash
# Download all files to a directory
# Make setup script executable
chmod +x setup_and_run.sh

# Install all dependencies
./setup_and_run.sh setup --with-nextflow
```

### Manual Setup

1. **Install Rust**:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env
```

2. **Install Python dependencies**:
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
python3 -m venv venv
source venv/bin/activate
pip install pandas numpy matplotlib seaborn scikit-learn scipy
```

3. **Install Nextflow** (optional):
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## 🏃‍♂️ Usage

### Method 1: Automated Script (Recommended)

```bash
# Quick pipeline run with default parameters
./setup_and_run.sh run-pipeline

# Custom parameters
./setup_and_run.sh run-pipeline --positions 10000 --states 4 --iterations 200 --parallel

# Individual steps
./setup_and_run.sh generate-data --positions 8000 --states 3
./setup_and_run.sh build
./setup_and_run.sh run-segmentation --states 3 --iterations 100
./setup_and_run.sh evaluate
```

### Method 2: Manual Execution

1. **Generate synthetic data**:
```bash
source venv/bin/activate
python3 generate_data.py --n-positions 5000 --n-states 3 --add-noise --visualize
```

2. **Build Rust tool**:
```bash
cargo build --release
```

3. **Run segmentation**:
```bash
./target/release/main \
    --input epigenomic_data.csv \
    --output segmentation_results.json \
    --states 3 \
    --max-iterations 100
```

4. **Evaluate results**:
```bash
python3 evaluate_segmentation.py \
    --results segmentation_results.json \
    --truth epigenomic_data_with_true_states.csv \
    --output-report evaluation_report.html
```

### Method 3: Nextflow Pipeline

```bash
# Basic run
nextflow run main.nf

# With custom parameters
nextflow run main.nf \
    --n_positions 8000 \
    --n_states 4 \
    --max_iterations 150 \
    --add_noise true \
    --visualize true \
    --parallel true
```

## 📊 Parameters

### Data Generation
- `--n-positions`: Number of genomic positions (default: 5000)
- `--n-states`: Number of chromatin states (default: 3)
- `--min-segment-length`: Minimum segment length (default: 50)
- `--noise-level`: Amount of noise to add (default: 0.1)
- `--add-noise`: Add realistic artifacts
- `--visualize`: Create visualization plots

### HMM Segmentation
- `--states`: Number of hidden states (default: 3)
- `--max-iterations`: Maximum EM iterations (default: 100)
- `--tolerance`: Convergence threshold (default: 1e-6)
- `--parallel`: Enable parallel processing

## 📈 Algorithm Details

### Hidden Markov Model
- **States**: Chromatin states (Heterochromatin, Euchromatin, Active Promoter, etc.)
- **Observations**: Continuous signal values (e.g., ChIP-seq, ATAC-seq)
- **Emissions**: Gaussian distributions with state-specific parameters
- **Transitions**: State-to-state transition probabilities

### Training Algorithm
1. **Initialization**: K-means-based parameter initialization
2. **E-step**: Forward-backward algorithm for posterior probabilities
3. **M-step**: Maximum likelihood parameter updates
4. **Convergence**: Log-likelihood improvement threshold

### State Decoding
- **Viterbi Algorithm**: Most likely state sequence
- **Posterior Decoding**: State probabilities at each position

## 🎯 Evaluation Metrics

### Clustering Metrics
- **Adjusted Rand Index (ARI)**: Clustering similarity (0-1, higher better)
- **Normalized Mutual Information (NMI)**: Information overlap (0-1, higher better)

### Classification Metrics
- **Accuracy**: Percentage of correctly classified positions
- **Precision/Recall/F1**: Per-state and macro-averaged metrics
- **Confusion Matrix**: Detailed classification breakdown

### Model Quality
- **Log-likelihood**: Model fit to data
- **Convergence**: Number of EM iterations
- **State Alignment**: Hungarian algorithm for label matching

## 📁 Output Files

### Generated Data
- `epigenomic_data.csv`: Input data for HMM (chromosome, position, signal)
- `epigenomic_data_with_true_states.csv`: Data with ground truth states
- `epigenomic_data_visualization.png`: Data visualization plots

### Segmentation Results
- `segmentation_results.json`: Complete HMM results including:
  - Predicted state sequence
  - Posterior probabilities
  - Model parameters
  - Training statistics

### Evaluation Outputs
- `evaluation_report.html`: Comprehensive evaluation report
- `evaluation_metrics.json`: Detailed metrics in JSON format
- `comparison_plots.png`: Visualization comparisons
- `confusion_matrix.png`: Confusion matrix heatmap
- `posterior_probs.png`: Posterior probability plots

### Pipeline Reports
- `pipeline_report.html`: Complete pipeline summary
- `pipeline_summary.json`: Machine-readable pipeline results

## 🔍 Example Results

### Typical Performance
- **Accuracy**: 85-95% on synthetic data
- **ARI**: 0.7-0.9 for well-separated states
- **Processing Time**: ~100-500ms for 5K positions
- **Convergence**: Usually 10-30 EM iterations

### State Examples
```
State 0 (Heterochromatin): Low signal, mean=0.5, var=0.3
State 1 (Euchromatin):     Medium signal, mean=2.0, var=0.5  
State 2 (Active Promoter): High signal, mean=4.5, var=0.8
```

## 🚀 Advanced Usage

### Large Datasets
```bash
# For large datasets (>50K positions)
./setup_and_run.sh run-pipeline \
    --positions 100000 \
    --states 5 \
    --parallel \
    --iterations 200
```

### Parameter Sensitivity
```bash
# Test different state numbers
for states in 2 3 4 5; do
    ./setup_and_run.sh run-pipeline --states $states --positions 5000
done
```

### Custom Data Format
The tool expects CSV format:
```csv
chromosome,position,signal
chr1,1000000,1.234
chr1,1001000,2.567
...
```

## 🐛 Troubleshooting

### Common Issues

1. **Rust Build Fails**:
```bash
# Update Rust
rustup update
# Clean and rebuild
cargo clean && cargo build --release
```

2. **Python Dependencies**:
```bash
# Reinstall in virtual environment
rm -rf venv
python3 -m venv venv
source venv/bin/activate
pip install pandas numpy matplotlib seaborn scikit-learn scipy
```

3. **Permission Denied**:
```bash
chmod +x setup_and_run.sh
chmod +x target/release/main
```

4. **WSL Issues**:
```bash
# Update WSL
wsl --update
# Install build essentials
sudo apt install build-essential
```

### Performance Issues
- Use `--parallel` flag for large datasets
- Reduce `--max-iterations` if convergence is slow
- Consider fewer states for faster processing

### Memory Issues
- Reduce `--n-positions` for large datasets
- Use streaming processing for very large files

## 📚 Technical References

### Algorithms
- **Forward-Backward Algorithm**: Baum-Welch for HMM training
- **Viterbi Algorithm**: Dynamic programming for state decoding
- **EM Algorithm**: Expectation-Maximization for parameter estimation
- **Hungarian Algorithm**: Optimal state alignment for evaluation

### Implementation Details
- **Language**: Rust 2021 edition
- **Dependencies**: serde, ndarray, statrs, rayon for parallelization
- **Numerical Stability**: Log-space computations, underflow protection
- **Memory Efficiency**: Streaming algorithms where possible

## 🤝 Contributing

### Development Setup
```bash
# Clone and setup development environment
git clone <repository>
cd epigenomic-hmm-segmentation
./setup_and_run.sh setup
```

### Code Style
- Follow Rust standard formatting: `cargo fmt`
- Run tests: `cargo test`
- Check with clippy: `cargo clippy`

### Adding Features
1. Extend the `StateParams` struct for new emission distributions
2. Modify the `forward_backward` function for new algorithms
3. Add evaluation metrics in `evaluate_segmentation.py`

## 📄 License

This project is released under the MIT License. See LICENSE file for details.

## 🙋‍♂️ Support

### Getting Help
1. Check this README for common solutions
2. Review log files in the results directory
3. Use `./setup_and_run.sh help` for command options

### Reporting Issues
When reporting issues, please include:
- Operating system and version (preferably WSL2)
- Error messages and log files
- Input parameters used
- Expected vs actual behavior

## 🔮 Future Enhancements

### Planned Features
- [ ] Multi-track segmentation
- [ ] Different emission distributions (Poisson, Negative Binomial)
- [ ] Hierarchical HMMs
- [ ] GPU acceleration
- [ ] Real genomic data integration
- [ ] Interactive visualization dashboard

### Performance Optimizations
- [ ] SIMD vectorization
- [ ] Sparse matrix operations
- [ ] Distributed computing support
- [ ] Streaming algorithms for large datasets

---

**🧬 Happy Segmenting!** 

For questions or contributions, please refer to the support section above.
