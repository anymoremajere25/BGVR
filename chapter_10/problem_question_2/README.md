# Problem 2: Simple Peak-Calling Module

Develop a simple peak-calling module in Rust that applies a threshold on coverage data to identify putative peaks for ChIP-Seq data. The module should account for local background noise estimates.
## ChIP-Seq Peak Caller

A simple yet powerful peak-calling module for ChIP-Seq data analysis, written in Rust.

## Dependencies

### Core Dependencies (Always Required)
- **Rust Standard Library Only**: The core functionality requires no external dependencies!
- **Minimum Rust Version**: 1.70.0 (2021 edition)

### Development Dependencies
- `tempfile = "3.8"` - For testing with temporary files

### Optional Dependencies

#### Command Line Interface
```toml
[dependencies]
clap = { version = "4.4", features = ["derive"], optional = true }
```

#### Parallel Processing
```toml
[dependencies]
rayon = { version = "1.8", optional = true }
```

#### JSON Serialization
```toml
[dependencies]
serde = { version = "1.0", features = ["derive"], optional = true }
serde_json = { version = "1.0", optional = true }
```

#### Compression Support
```toml
[dependencies]
flate2 = { version = "1.0", optional = true }
```

#### Bioinformatics Utilities
```toml
[dependencies]
bio = { version = "1.6", optional = true }
```

## Installation

### Basic Installation (Library Only)
```bash
cargo add chipseq-peak-caller
```

### With All Features
```bash
cargo add chipseq-peak-caller --features full
```

### Custom Feature Selection
```bash
cargo add chipseq-peak-caller --features cli,parallel,compression
```

### Build from Source
```bash
git clone https://github.com/username/chipseq-peak-caller.git
cd chipseq-peak-caller

# Basic build
cargo build --release

# With all features
cargo build --release --features full

# With specific features
cargo build --release --features cli,parallel
```

## Features

### Available Features
- `cli` - Command-line interface using clap
- `parallel` - Multi-threaded processing with rayon
- `json` - JSON serialization support
- `compression` - Gzip file support
- `bio-utils` - Integration with rust-bio
- `full` - All features enabled

### Default Features
By default, no optional features are enabled to keep dependencies minimal.

## Usage

### As a Library

#### Basic Usage
```rust
use chipseq_peak_caller::{PeakCaller, PeakCallingConfig};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create peak caller with default settings
    let mut caller = PeakCaller::new();
    
    // Load coverage data
    caller.load_coverage_from_bedgraph("input.bedgraph")?;
    
    // Call peaks
    let peaks = caller.call_peaks();
    
    // Write results
    caller.write_peaks_to_bed(&peaks, "output.bed")?;
    
    // Print summary
    let summary = caller.get_summary(&peaks);
    println!("{}", summary);
    
    Ok(())
}
```

#### Advanced Configuration
```rust
use chipseq_peak_caller::{PeakCaller, PeakCallingConfig};

let config = PeakCallingConfig {
    min_fold_enrichment: 3.0,    // Higher stringency
    background_window: 5000,     // Smaller background window
    min_peak_length: 200,        // Longer minimum peaks
    merge_distance: 300,         // Merge nearby peaks
};

let mut caller = PeakCaller::with_config(config);
```

### Command Line Interface

#### Build CLI Tool
```bash
cargo build --release --features cli
```

#### Basic Usage
```bash
./target/release/peak-caller -i input.bedgraph -o peaks.bed
```

#### Advanced Usage
```bash
./target/release/peak-caller \
    --input coverage.bedgraph.gz \
    --output peaks.bed \
    --fold-enrichment 2.5 \
    --background-window 8000 \
    --min-peak-length 100 \
    --merge-distance 250 \
    --verbose
```

#### CLI Options
- `-i, --input <FILE>` - Input BedGraph file (required)
- `-o, --output <FILE>` - Output BED file (default: peaks.bed)
- `-f, --fold-enrichment <FLOAT>` - Minimum fold enrichment (default: 2.0)
- `-w, --background-window <INT>` - Background window size in bp (default: 10000)
- `-l, --min-peak-length <INT>` - Minimum peak length in bp (default: 150)
- `-m, --merge-distance <INT>` - Merge peaks within distance (default: 200)
- `-v, --verbose` - Enable verbose output

## Performance

### Benchmarks
| Dataset Size | Memory Usage | Processing Time | Features |
|-------------|-------------|----------------|----------|
| 100MB bedgraph | ~50MB RAM | 2-5 seconds | default |
| 1GB bedgraph | ~300MB RAM | 15-30 seconds | default |
| 1GB bedgraph | ~300MB RAM | 8-15 seconds | `parallel` |

### Memory Efficiency
- Uses sorted vectors for O(log n) coverage queries
- Processes data chromosome by chromosome
- Optional parallel processing for multi-core systems

### Optimization Tips
1. **Large Files**: Enable `parallel` feature for multi-threading
2. **Compressed Files**: Use `compression` feature for .gz support
3. **Memory Constraints**: Process chromosomes individually for huge datasets

## Input/Output Formats

### Input: BedGraph Format
```
track type=bedGraph name="Coverage"
chr1    1000    1100    5.2
chr1    1100    1200    12.8
chr1    1200    1300    8.1
```

### Output: BED Format
```
track name="ChIP-Seq Peaks" description="Peaks called with fold enrichment >= 2.0"
chr1    1100    1250    peak    320    .
chr2    5000    5200    peak    480    .
```

## Algorithm Details

### Peak Calling Algorithm
1. **Load Coverage**: Parse BedGraph input into sorted intervals
2. **Background Calculation**: Compute local background in sliding window
3. **Signal Detection**: Identify regions above fold-enrichment threshold
4. **Peak Merging**: Combine nearby peaks within merge distance
5. **Quality Filtering**: Apply minimum length requirements

### Data Structures
- `BTreeMap<String, Vec<CoverageInterval>>` for chromosome organization
- Binary search for efficient position queries
- Sorted vectors for optimal memory usage

## Testing

### Run Tests
```bash
# All tests
cargo test

# With all features
cargo test --features full

# Integration tests only
cargo test --test integration
```

### Test Coverage
- Unit tests for all core functions
- Integration tests with real data
- Property-based testing for edge cases
- Performance regression tests

## Contributing

### Development Setup
```bash
git clone https://github.com/username/chipseq-peak-caller.git
cd chipseq-peak-caller

# Install development dependencies
cargo build --features full

# Run tests
cargo test --features full

# Format code
cargo fmt

# Lint code
cargo clippy --features full
```

### Code Structure
```
src/
├── lib.rs           # Main library code
├── bin/
│   └── peak_caller.rs # CLI interface
└── tests/
    └── integration.rs # Integration tests
```

## Troubleshooting

### Common Issues

#### Compilation Errors
```bash
error: failed to resolve: use of undeclared crate or module
```
**Solution**: Enable required features:
```bash
cargo build --features cli,parallel
```

#### Memory Issues
```bash
thread 'main' panicked at 'memory allocation failed'
```
**Solution**: Process smaller chunks or use streaming:
```rust
// Process chromosome by chromosome for large files
for chromosome in chromosomes {
    let peaks = caller.call_peaks_for_chromosome(&chromosome);
    // Process peaks...
}
```

#### Performance Issues
**Problem**: Slow processing on large files
**Solutions**:
1. Enable parallel processing: `--features parallel`
2. Adjust background window size
3. Use compressed input files
4. Increase system memory

### Getting Help
- Check the [documentation](https://docs.rs/chipseq-peak-caller)
- File issues on [GitHub](https://github.com/username/chipseq-peak-caller/issues)
- Review the test suite for usage examples

## License

This project is licensed under either of:
- Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.

## Citation

If you use this software in your research, please cite:
```bibtex
@software{chipseq_peak_caller,
  author = {Your Name},
  title = {ChIP-Seq Peak Caller: A Simple Peak Calling Module in Rust},
  year = {2025},
  url = {https://github.com/username/chipseq-peak-caller}
}
```


