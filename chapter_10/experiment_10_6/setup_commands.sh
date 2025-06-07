#!/bin/bash

# Epigenomic HMM Segmentation - Complete Setup and Execution Script
# This script sets up the environment and runs the entire pipeline

set -e  # Exit on any error

echo "🧬 Epigenomic HMM Segmentation Pipeline Setup"
echo "============================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

print_step() {
    echo -e "\n${BLUE}==== $1 ====${NC}"
}

# Check if we're in WSL
check_wsl() {
    if grep -qEi "(Microsoft|WSL)" /proc/version &> /dev/null; then
        print_status "Running in WSL environment"
        return 0
    else
        print_warning "Not running in WSL - this script is optimized for WSL"
        return 1
    fi
}

# Install dependencies
install_dependencies() {
    print_step "Installing Dependencies"
    
    # Update package list
    print_status "Updating package list..."
    sudo apt update
    
    # Install Rust if not present
    if ! command -v rustc &> /dev/null; then
        print_status "Installing Rust..."
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
        source ~/.cargo/env
    else
        print_status "Rust already installed: $(rustc --version)"
    fi
    
    # Install Python dependencies
    print_status "Installing Python dependencies..."
    sudo apt install -y python3 python3-pip python3-venv
    
    # Create virtual environment
    if [ ! -d "venv" ]; then
        print_status "Creating Python virtual environment..."
        python3 -m venv venv
    fi
    
    # Activate virtual environment and install packages
    print_status "Installing Python packages..."
    source venv/bin/activate
    pip install --upgrade pip
    pip install pandas numpy matplotlib seaborn scikit-learn scipy
    
    # Install Nextflow if requested
    if [ "$1" = "--with-nextflow" ]; then
        install_nextflow
    fi
}

install_nextflow() {
    print_status "Installing Nextflow..."
    if ! command -v nextflow &> /dev/null; then
        curl -s https://get.nextflow.io | bash
        sudo mv nextflow /usr/local/bin/
        print_status "Nextflow installed successfully"
    else
        print_status "Nextflow already installed: $(nextflow -version | head -1)"
    fi
}

# Setup project structure
setup_project() {
    print_step "Setting Up Project Structure"
    
    # Create directories
    mkdir -p src data results
    
    # Move files to appropriate locations
    if [ -f "main.rs" ]; then
        mv main.rs src/
        print_status "Moved main.rs to src/"
    fi
    
    print_status "Project structure created"
    ls -la
}

# Build Rust project
build_rust() {
    print_step "Building Rust Project"
    
    # Initialize Cargo project if needed
    if [ ! -f "Cargo.toml" ]; then
        print_error "Cargo.toml not found! Make sure you have the Cargo.toml file."
        exit 1
    fi
    
    # Build the project
    print_status "Building Rust project..."
    cargo build --release
    
    if [ -f "target/release/main" ]; then
        print_status "✅ Rust build successful!"
    else
        print_error "❌ Rust build failed!"
        exit 1
    fi
}

# Generate synthetic data
generate_data() {
    print_step "Generating Synthetic Data"
    
    local n_positions=${1:-5000}
    local n_states=${2:-3}
    local add_noise=${3:-"--add-noise"}
    local visualize=${4:-"--visualize"}
    
    # Activate virtual environment
    source venv/bin/activate
    
    print_status "Generating $n_positions positions with $n_states states..."
    python3 generate_data.py \
        --n-positions $n_positions \
        --n-states $n_states \
        --output data/epigenomic_data.csv \
        $add_noise \
        $visualize
    
    if [ -f "data/epigenomic_data.csv" ]; then
        print_status "✅ Data generation successful!"
        print_status "Files created:"
        ls -la data/
    else
        print_error "❌ Data generation failed!"
        exit 1
    fi
}

# Run HMM segmentation
run_segmentation() {
    print_step "Running HMM Segmentation"
    
    local n_states=${1:-3}
    local max_iterations=${2:-100}
    local parallel=${3:-""}
    
    print_status "Running segmentation with $n_states states, max $max_iterations iterations..."
    
    ./target/release/main \
        --input data/epigenomic_data.csv \
        --output results/segmentation_results.json \
        --states $n_states \
        --max-iterations $max_iterations \
        $parallel
    
    if [ -f "results/segmentation_results.json" ]; then
        print_status "✅ HMM segmentation successful!"
    else
        print_error "❌ HMM segmentation failed!"
        exit 1
    fi
}

# Evaluate results
evaluate_results() {
    print_step "Evaluating Results"
    
    # Activate virtual environment
    source venv/bin/activate
    
    print_status "Evaluating segmentation results..."
    python3 evaluate_segmentation.py \
        --results results/segmentation_results.json \
        --truth data/epigenomic_data_with_true_states.csv \
        --output-report results/evaluation_report.html \
        --output-metrics results/evaluation_metrics.json \
        --output-plots results/comparison_plots.png
    
    if [ -f "results/evaluation_report.html" ]; then
        print_status "✅ Evaluation successful!"
    else
        print_error "❌ Evaluation failed!"
        exit 1
    fi
}

# Run full pipeline
run_full_pipeline() {
    print_step "Running Complete Pipeline"
    
    local n_positions=${1:-5000}
    local n_states=${2:-3}
    local max_iterations=${3:-100}
    local parallel_flag=""
    
    if [ "$4" = "--parallel" ]; then
        parallel_flag="--parallel"
    fi
    
    print_status "Pipeline parameters:"
    print_status "  Positions: $n_positions"
    print_status "  States: $n_states"
    print_status "  Max iterations: $max_iterations"
    print_status "  Parallel: ${parallel_flag:-"disabled"}"
    
    # Run all steps
    generate_data $n_positions $n_states "--add-noise" "--visualize"
    run_segmentation $n_states $max_iterations "$parallel_flag"
    evaluate_results
    
    print_status "✅ Full pipeline completed successfully!"
    print_status "Results available in: results/"
    ls -la results/
}

# Run with Nextflow
run_nextflow() {
    print_step "Running Nextflow Pipeline"
    
    local n_positions=${1:-5000}
    local n_states=${2:-3}
    local max_iterations=${3:-100}
    
    if ! command -v nextflow &> /dev/null; then
        print_error "Nextflow not installed. Run with --install-nextflow first."
        exit 1
    fi
    
    print_status "Running Nextflow pipeline..."
    nextflow run main.nf \
        --n_positions $n_positions \
        --n_states $n_states \
        --max_iterations $max_iterations \
        --add_noise true \
        --visualize true \
        --outdir results
    
    print_status "✅ Nextflow pipeline completed!"
}

# Clean up
cleanup() {
    print_step "Cleaning Up"
    
    print_status "Removing build artifacts..."
    cargo clean 2>/dev/null || true
    
    print_status "Removing temporary files..."
    rm -f *.log
    
    print_status "✅ Cleanup completed!"
}

# Show results
show_results() {
    print_step "Results Summary"
    
    if [ -f "results/segmentation_results.json" ]; then
        print_status "Segmentation results:"
        echo "  📄 results/segmentation_results.json"
        
        # Extract basic info from JSON
        if command -v python3 &> /dev/null; then
            python3 -c "
import json
try:
    with open('results/segmentation_results.json', 'r') as f:
        data = json.load(f)
    print(f'  ⏱️  Processing time: {data.get(\"processing_time_ms\", \"N/A\")} ms')
    print(f'  🔄 Model iterations: {data.get(\"model\", {}).get(\"iterations\", \"N/A\")}')
    print(f'  📊 Data length: {data.get(\"data_length\", \"N/A\")}')
    print(f'  🎯 Number of states: {data.get(\"model\", {}).get(\"num_states\", \"N/A\")}')
except Exception as e:
    print(f'  ❌ Error reading results: {e}')
"
        fi
    fi
    
    if [ -f "results/evaluation_metrics.json" ]; then
        print_status "Evaluation metrics:"
        echo "  📄 results/evaluation_metrics.json"
        echo "  📄 results/evaluation_report.html"
        
        if command -v python3 &> /dev/null; then
            python3 -c "
import json
try:
    with open('results/evaluation_metrics.json', 'r') as f:
        data = json.load(f)
    print(f'  🎯 Accuracy: {data.get(\"accuracy\", \"N/A\"):.3f}' if isinstance(data.get('accuracy'), (int, float)) else '  🎯 Accuracy: N/A')
    print(f'  📈 Adjusted Rand Score: {data.get(\"adjusted_rand_score\", \"N/A\"):.3f}' if isinstance(data.get('adjusted_rand_score'), (int, float)) else '  📈 Adjusted Rand Score: N/A')
    print(f'  🔍 F1 Score: {data.get(\"f1_macro\", \"N/A\"):.3f}' if isinstance(data.get('f1_macro'), (int, float)) else '  🔍 F1 Score: N/A')
except Exception as e:
    print(f'  ❌ Error reading metrics: {e}')
"
        fi
    fi
    
    print_status "Generated files:"
    find data results -type f 2>/dev/null | sort | while read file; do
        echo "  📄 $file"
    done
    
    if [ -f "results/evaluation_report.html" ]; then
        print_status "🌐 Open results/evaluation_report.html in a web browser to view detailed results"
    fi
}

# Help function
show_help() {
    echo "Epigenomic HMM Segmentation Pipeline"
    echo "Usage: $0 [COMMAND] [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  setup              - Install dependencies and setup project"
    echo "  build              - Build Rust project only"
    echo "  generate-data      - Generate synthetic data only"
    echo "  run-segmentation   - Run HMM segmentation only"
    echo "  evaluate           - Evaluate results only"
    echo "  run-pipeline       - Run complete pipeline"
    echo "  run-nextflow       - Run Nextflow pipeline"
    echo "  cleanup            - Clean up build artifacts"
    echo "  show-results       - Show results summary"
    echo "  help               - Show this help"
    echo ""
    echo "Options for setup:"
    echo "  --with-nextflow    - Also install Nextflow"
    echo ""
    echo "Options for generate-data:"
    echo "  --positions N      - Number of positions (default: 5000)"
    echo "  --states N         - Number of states (default: 3)"
    echo "  --no-noise         - Don't add realistic noise"
    echo "  --no-visualize     - Don't create visualizations"
    echo ""
    echo "Options for run-segmentation:"
    echo "  --states N         - Number of states (default: 3)"
    echo "  --iterations N     - Max iterations (default: 100)"
    echo "  --parallel         - Use parallel processing"
    echo ""
    echo "Examples:"
    echo "  $0 setup --with-nextflow"
    echo "  $0 generate-data --positions 10000 --states 4"
    echo "  $0 run-pipeline --positions 5000 --states 3 --iterations 50"
    echo "  $0 run-nextflow --positions 8000 --states 4"
    echo ""
    echo "Quick start:"
    echo "  $0 setup && $0 run-pipeline"
}

# Parse command line arguments
COMMAND=${1:-help}
shift || true

case $COMMAND in
    setup)
        check_wsl
        install_dependencies "$@"
        setup_project
        ;;
    build)
        build_rust
        ;;
    generate-data)
        # Parse arguments
        n_positions=5000
        n_states=3
        add_noise="--add-noise"
        visualize="--visualize"
        
        while [[ $# -gt 0 ]]; do
            case $1 in
                --positions)
                    n_positions="$2"
                    shift 2
                    ;;
                --states)
                    n_states="$2"
                    shift 2
                    ;;
                --no-noise)
                    add_noise=""
                    shift
                    ;;
                --no-visualize)
                    visualize=""
                    shift
                    ;;
                *)
                    shift
                    ;;
            esac
        done
        
        generate_data $n_positions $n_states "$add_noise" "$visualize"
        ;;
    run-segmentation)
        # Parse arguments
        n_states=3
        max_iterations=100
        parallel=""
        
        while [[ $# -gt 0 ]]; do
            case $1 in
                --states)
                    n_states="$2"
                    shift 2
                    ;;
                --iterations)
                    max_iterations="$2"
                    shift 2
                    ;;
                --parallel)
                    parallel="--parallel"
                    shift
                    ;;
                *)
                    shift
                    ;;
            esac
        done
        
        run_segmentation $n_states $max_iterations "$parallel"
        ;;
    evaluate)
        evaluate_results
        ;;
    run-pipeline)
        # Parse arguments
        n_positions=5000
        n_states=3
        max_iterations=100
        parallel=""
        
        while [[ $# -gt 0 ]]; do
            case $1 in
                --positions)
                    n_positions="$2"
                    shift 2
                    ;;
                --states)
                    n_states="$2"
                    shift 2
                    ;;
                --iterations)
                    max_iterations="$2"
                    shift 2
                    ;;
                --parallel)
                    parallel="--parallel"
                    shift
                    ;;
                *)
                    shift
                    ;;
            esac
        done
        
        build_rust
        run_full_pipeline $n_positions $n_states $max_iterations "$parallel"
        show_results
        ;;
    run-nextflow)
        # Parse arguments
        n_positions=5000
        n_states=3
        max_iterations=100
        
        while [[ $# -gt 0 ]]; do
            case $1 in
                --positions)
                    n_positions="$2"
                    shift 2
                    ;;
                --states)
                    n_states="$2"
                    shift 2
                    ;;
                --iterations)
                    max_iterations="$2"
                    shift 2
                    ;;
                *)
                    shift
                    ;;
            esac
        done
        
        run_nextflow $n_positions $n_states $max_iterations
        ;;
    cleanup)
        cleanup
        ;;
    show-results)
        show_results
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        print_error "Unknown command: $COMMAND"
        echo ""
        show_help
        exit 1
        ;;
esac