#!/usr/bin/env python3
"""
Enhanced Synthetic Dataset Generator for Interval-based Peak Calling
Generates realistic test data with predefined genomic intervals and enriched regions
"""

import os
import gzip
import random
import argparse
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime

class EnhancedGenomeGenerator:
    """Generate synthetic reference genome with realistic features"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
        
    def generate_realistic_sequence(self, length: int, gc_content: float = 0.4, 
                                  add_repeats: bool = True) -> str:
        """Generate more realistic genomic sequence with repeats and structure"""
        sequence = []
        
        # Generate basic sequence
        for i in range(length):
            # Create GC islands occasionally
            if i % 2000 < 200 and random.random() < 0.3:  # GC island
                local_gc = min(0.7, gc_content + 0.3)
            else:
                local_gc = gc_content
                
            if random.random() < local_gc:
                sequence.append(random.choice(['G', 'C']))
            else:
                sequence.append(random.choice(['A', 'T']))
        
        if add_repeats:
            sequence = self.add_repeat_elements(sequence)
            
        return ''.join(sequence)
    
    def add_repeat_elements(self, sequence: List[str]) -> List[str]:
        """Add simple repeat elements to make sequence more realistic"""
        seq_len = len(sequence)
        
        # Add some tandem repeats
        for _ in range(seq_len // 10000):  # About every 10kb
            if seq_len < 100:
                break
                
            start = random.randint(0, seq_len - 100)
            repeat_unit = ''.join(random.choices(['A', 'T', 'G', 'C'], k=random.randint(2, 6)))
            repeat_length = random.randint(10, 50)
            
            # Replace region with tandem repeat
            end = min(start + repeat_length, seq_len)
            for i in range(start, end):
                sequence[i] = repeat_unit[i % len(repeat_unit)]
        
        return sequence
    
    def generate_genome(self, chromosomes: List[Tuple[str, int]], output_file: str):
        """Generate multi-chromosome reference genome with realistic features"""
        print(f"Generating enhanced reference genome with {len(chromosomes)} chromosomes...")
        
        with open(output_file, 'w') as f:
            for chrom_name, length in chromosomes:
                print(f"  Generating {chrom_name} ({length:,} bp) with realistic features...")
                
                # Vary GC content by chromosome
                base_gc = 0.41 if 'chr' in chrom_name.lower() else 0.38
                gc_content = base_gc + random.uniform(-0.05, 0.05)
                
                f.write(f">{chrom_name}\n")
                
                # Generate sequence in chunks for memory efficiency
                chunk_size = 50000
                for start in range(0, length, chunk_size):
                    chunk_length = min(chunk_size, length - start)
                    chunk_seq = self.generate_realistic_sequence(chunk_length, gc_content)
                    
                    # Write in 80-character lines
                    for i in range(0, len(chunk_seq), 80):
                        f.write(chunk_seq[i:i+80] + "\n")
        
        print(f"Enhanced reference genome saved to {output_file}")

class IntervalGenerator:
    """Generate biologically relevant genomic intervals"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        
    def generate_gene_like_intervals(self, chromosomes: List[Tuple[str, int]], 
                                   genes_per_mb: float = 20) -> List[Tuple[str, int, int, str, str]]:
        """Generate gene-like intervals with realistic properties"""
        intervals = []
        
        gene_types = ['protein_coding', 'lncRNA', 'pseudogene', 'miRNA', 'enhancer']
        gene_type_weights = [0.6, 0.15, 0.15, 0.05, 0.05]
        
        for chrom_name, chrom_length in chromosomes:
            num_genes = int((chrom_length / 1_000_000) * genes_per_mb)
            
            # Generate non-overlapping gene positions
            gene_positions = []
            for _ in range(num_genes * 2):  # Try more times to avoid overlaps
                gene_length = random.randint(1000, 50000)  # 1kb to 50kb genes
                start = random.randint(1000, chrom_length - gene_length - 1000)
                end = start + gene_length
                
                # Check for overlaps
                overlap = False
                for existing_start, existing_end in gene_positions:
                    if start < existing_end and end > existing_start:
                        overlap = True
                        break
                
                if not overlap:
                    gene_type = random.choices(gene_types, weights=gene_type_weights)[0]
                    gene_name = f"{gene_type}_{len(gene_positions) + 1}"
                    intervals.append((chrom_name, start, end, gene_name, gene_type))
                    gene_positions.append((start, end))
                    
                    if len(gene_positions) >= num_genes:
                        break
        
        print(f"Generated {len(intervals)} gene-like intervals")
        return intervals
    
    def generate_regulatory_intervals(self, gene_intervals: List[Tuple[str, int, int, str, str]],
                                    promoter_size: int = 2000,
                                    enhancer_density: float = 0.3) -> List[Tuple[str, int, int, str, str]]:
        """Generate promoter and enhancer intervals based on gene positions"""
        regulatory_intervals = []
        
        for chrom, gene_start, gene_end, gene_name, gene_type in gene_intervals:
            # Add promoter region
            promoter_start = max(0, gene_start - promoter_size)
            promoter_end = gene_start
            regulatory_intervals.append((
                chrom, promoter_start, promoter_end, 
                f"promoter_{gene_name}", "promoter"
            ))
            
            # Randomly add enhancers
            if random.random() < enhancer_density:
                # Enhancer can be upstream or downstream
                if random.random() < 0.5:  # Upstream
                    enhancer_start = max(0, gene_start - random.randint(5000, 20000))
                    enhancer_end = enhancer_start + random.randint(500, 2000)
                else:  # Downstream
                    enhancer_start = gene_end + random.randint(1000, 10000)
                    enhancer_end = enhancer_start + random.randint(500, 2000)
                
                regulatory_intervals.append((
                    chrom, enhancer_start, enhancer_end,
                    f"enhancer_{gene_name}_{random.randint(1, 10)}", "enhancer"
                ))
        
        print(f"Generated {len(regulatory_intervals)} regulatory intervals")
        return regulatory_intervals

class EnhancedFASTQGenerator:
    """Generate synthetic FASTQ reads with realistic ChIP-seq characteristics"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
        
    def generate_chip_seq_reads(self, genome_file: str, 
                               target_intervals: List[Tuple[str, int, int, str, str]],
                               num_reads: int, read_length: int = 75,
                               enrichment_factor: float = 10.0,
                               fragment_size_mean: int = 200,
                               fragment_size_std: int = 50) -> List[Tuple[str, str, str]]:
        """Generate ChIP-seq-like reads with fragment size distribution and enrichment"""
        print(f"Generating {num_reads:,} ChIP-seq reads...")
        
        # Load genome sequences
        genome_seqs = self.load_genome_sequences(genome_file)
        
        reads = []
        target_reads = 0
        background_reads = 0
        
        # Create enrichment map
        enrichment_regions = []
        for chrom, start, end, name, region_type in target_intervals:
            if region_type in ['promoter', 'enhancer']:  # Only enrich certain types
                # Add some noise to exact coordinates
                noise_start = max(0, start - random.randint(0, 200))
                noise_end = end + random.randint(0, 200)
                enrichment_regions.append((chrom, noise_start, noise_end, enrichment_factor))
        
        print(f"Created {len(enrichment_regions)} enrichment regions")
        
        for read_id in range(num_reads):
            # Decide if this read comes from an enriched region
            if random.random() < 0.4 and enrichment_regions:  # 40% from enriched regions
                # Sample from enriched regions
                chrom, start, end, enrichment = random.choice(enrichment_regions)
                
                if chrom in genome_seqs:
                    # Generate realistic fragment
                    fragment_size = max(read_length + 20, 
                                      int(np.random.normal(fragment_size_mean, fragment_size_std)))
                    
                    # Sample fragment position
                    region_length = end - start
                    if region_length > fragment_size:
                        fragment_start = start + random.randint(0, region_length - fragment_size)
                        fragment_end = fragment_start + fragment_size
                        
                        # Choose read position within fragment (5' bias)
                        if random.random() < 0.7:  # 5' end bias
                            read_start = fragment_start + random.randint(0, min(50, fragment_size - read_length))
                        else:  # 3' end
                            read_start = fragment_end - read_length - random.randint(0, 50)
                        
                        read_start = max(0, min(read_start, len(genome_seqs[chrom]) - read_length))
                        
                        if read_start + read_length <= len(genome_seqs[chrom]):
                            read_seq = genome_seqs[chrom][read_start:read_start + read_length]
                            read_seq = self.add_sequencing_errors(read_seq, error_rate=0.008)
                            quality = self.generate_realistic_quality(read_length, base_quality=32)
                            
                            reads.append((f"read_{read_id}_enriched", read_seq, quality))
                            target_reads += 1
                            continue
            
            # Background read
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            
            if len(chrom_seq) > read_length:
                read_start = random.randint(0, len(chrom_seq) - read_length)
                read_seq = chrom_seq[read_start:read_start + read_length]
                read_seq = self.add_sequencing_errors(read_seq, error_rate=0.01)
                quality = self.generate_realistic_quality(read_length, base_quality=28)
                
                reads.append((f"read_{read_id}_background", read_seq, quality))
                background_reads += 1
        
        print(f"Generated {target_reads:,} enriched reads and {background_reads:,} background reads")
        return reads
    
    def load_genome_sequences(self, genome_file: str) -> Dict[str, str]:
        """Load genome sequences from FASTA file"""
        genome_seqs = {}
        current_chrom = None
        current_seq = []
        
        with open(genome_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_chrom:
                        genome_seqs[current_chrom] = ''.join(current_seq)
                    current_chrom = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_chrom:
                genome_seqs[current_chrom] = ''.join(current_seq)
        
        print(f"Loaded {len(genome_seqs)} chromosomes from genome")
        return genome_seqs
    
    def add_sequencing_errors(self, sequence: str, error_rate: float = 0.01) -> str:
        """Add realistic sequencing errors"""
        mutated = list(sequence)
        for i in range(len(mutated)):
            if random.random() < error_rate:
                error_type = random.choice(['substitution', 'deletion', 'insertion'])
                if error_type == 'substitution':
                    mutated[i] = random.choice([b for b in self.bases if b != mutated[i]])
                elif error_type == 'deletion' and len(mutated) > 1:
                    mutated.pop(i)
                    break
                elif error_type == 'insertion':
                    mutated.insert(i, random.choice(self.bases))
                    break
        return ''.join(mutated)
    
    def generate_realistic_quality(self, length: int, base_quality: int = 30) -> str:
        """Generate realistic quality scores with position-dependent degradation"""
        qualities = []
        for i in range(length):
            # Quality decreases towards end with some randomness
            pos_factor = 1.0 - (i / length) * 0.4
            cycle_noise = 0.1 * np.sin(i * 0.1)  # Cycle-dependent noise
            random_noise = random.gauss(0, 2)
            
            quality = base_quality * pos_factor + cycle_noise + random_noise
            quality = max(15, min(40, int(quality)))  # Clamp between Q15 and Q40
            qualities.append(chr(quality + 33))
        return ''.join(qualities)
    
    def write_fastq(self, reads: List[Tuple[str, str, str]], output_file: str):
        """Write reads to compressed FASTQ file"""
        print(f"Writing {len(reads):,} reads to {output_file}")
        
        with gzip.open(output_file, 'wt') as f:
            for read_id, sequence, quality in reads:
                # Ensure sequence and quality have same length
                min_len = min(len(sequence), len(quality))
                f.write(f"@{read_id}\n")
                f.write(f"{sequence[:min_len]}\n")
                f.write(f"+\n")
                f.write(f"{quality[:min_len]}\n")

def create_directory_structure(base_dir: str):
    """Create comprehensive directory structure"""
    directories = [
        f"{base_dir}/data",
        f"{base_dir}/reference",
        f"{base_dir}/intervals",
        f"{base_dir}/results",
        f"{base_dir}/logs",
        f"{base_dir}/config"
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
    
    print(f"Created directory structure in {base_dir}")

def write_intervals_bed(intervals: List[Tuple[str, int, int, str, str]], output_file: str):
    """Write intervals to BED format"""
    with open(output_file, 'w') as f:
        f.write("# chrom\tstart\tend\tname\ttype\n")
        for chrom, start, end, name, interval_type in intervals:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{interval_type}\n")
    
    print(f"Wrote {len(intervals)} intervals to {output_file}")

def write_true_peaks_bed(intervals: List[Tuple[str, int, int, str, str]], output_file: str):
    """Write enriched regions as true peaks for validation"""
    enriched_regions = [
        (chrom, start, end, name, interval_type) 
        for chrom, start, end, name, interval_type in intervals
        if interval_type in ['promoter', 'enhancer']
    ]
    
    with open(output_file, 'w') as f:
        f.write("# chrom\tstart\tend\tname\ttype\tenrichment\n")
        for chrom, start, end, name, interval_type in enriched_regions:
            enrichment = random.uniform(5.0, 20.0)  # Realistic enrichment values
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{interval_type}\t{enrichment:.2f}\n")
    
    print(f"Wrote {len(enriched_regions)} true enriched regions to {output_file}")

def generate_nextflow_config(base_dir: str, chromosomes: List[Tuple[str, int]]):
    """Generate enhanced Nextflow configuration"""
    total_genome_size = sum(length for _, length in chromosomes)
    
    config_content = f"""// Enhanced Interval-based Pipeline Configuration
// Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

params {{
    // Input/Output directories
    input_dir = "{base_dir}/data"
    output_dir = "{base_dir}/results"
    reference_genome = "{base_dir}/reference/genome.fa"
    intervals_bed = "{base_dir}/intervals/gene_intervals.bed"
    
    // Analysis parameters
    window_size = 500
    min_coverage = 3.0
    pvalue_threshold = 0.05
    fragment_shift = 75
    extend_reads = 150
    
    // Quality control
    min_quality = 20
    mapping_quality = 10
    remove_duplicates = false
    
    // Performance
    threads = 4
    cache_results = true
    output_bed = true
    
    // Advanced options
    paired_end = false
    strand_specific = false
}}

// Resource configuration
process {{
    // Default resources
    cpus = 2
    memory = '4 GB'
    time = '2h'
    
    withName: BUILD_INDICES {{
        cpus = 4
        memory = '8 GB'
        time = '1h'
    }}
    
    withName: ALIGN_READS {{
        cpus = 4
        memory = '6 GB'
        time = '2h'
    }}
    
    withName: INTERVAL_PEAK_CALLING {{
        cpus = 4
        memory = '8 GB'
        time = '3h'
    }}
}}

// Execution configuration
executor {{
    name = 'local'
    cpus = 8
    memory = '16 GB'
}}

// Reporting
report {{
    enabled = true
    file = "{base_dir}/results/nextflow_report.html"
}}

timeline {{
    enabled = true
    file = "{base_dir}/results/nextflow_timeline.html"
}}

trace {{
    enabled = true
    file = "{base_dir}/results/nextflow_trace.txt"
}}

// Docker/Singularity (optional)
// docker {{
//     enabled = true
//     runOptions = '-u $(id -u):$(id -g)'
// }}

// Genome information (for reference)
// Total genome size: {total_genome_size:,} bp
// Number of chromosomes: {len(chromosomes)}
"""
    
    config_file = f"{base_dir}/nextflow.config"
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    print(f"Generated Nextflow configuration: {config_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate enhanced synthetic dataset for interval-based peak calling")
    parser.add_argument("--output-dir", default="./interval_test", 
                       help="Output directory for generated data")
    parser.add_argument("--num-reads", type=int, default=300000,
                       help="Number of reads per sample")
    parser.add_argument("--num-samples", type=int, default=3,
                       help="Number of sample FASTQ files")
    parser.add_argument("--read-length", type=int, default=75,
                       help="Length of generated reads")
    parser.add_argument("--enrichment-factor", type=float, default=8.0,
                       help="Enrichment factor for target regions")
    parser.add_argument("--genes-per-mb", type=float, default=25.0,
                       help="Number of gene-like intervals per megabase")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducible data")
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Enhanced Interval-based Dataset Generator")
    print("=" * 60)
    print(f"Output directory: {args.output_dir}")
    print(f"Reads per sample: {args.num_reads:,}")
    print(f"Number of samples: {args.num_samples}")
    print(f"Read length: {args.read_length} bp")
    print(f"Enrichment factor: {args.enrichment_factor}x")
    print(f"Genes per Mb: {args.genes_per_mb}")
    print(f"Random seed: {args.seed}")
    print()
    
    # Create directory structure
    create_directory_structure(args.output_dir)
    
    # Define realistic chromosome sizes for testing
    chromosomes = [
        ("chr1", 3_000_000),   # 3 Mb
        ("chr2", 2_500_000),   # 2.5 Mb  
        ("chr3", 2_000_000),   # 2 Mb
        ("chr4", 1_500_000),   # 1.5 Mb
        ("chrX", 1_000_000),   # 1 Mb
    ]
    
    # Generate enhanced reference genome
    print("\n--- Generating Enhanced Reference Genome ---")
    genome_gen = EnhancedGenomeGenerator(seed=args.seed)
    genome_file = f"{args.output_dir}/reference/genome.fa"
    genome_gen.generate_genome(chromosomes, genome_file)
    
    # Generate gene-like intervals
    print("\n--- Generating Gene-like Intervals ---")
    interval_gen = IntervalGenerator(seed=args.seed)
    gene_intervals = interval_gen.generate_gene_like_intervals(chromosomes, args.genes_per_mb)
    
    # Generate regulatory intervals
    print("\n--- Generating Regulatory Intervals ---")
    regulatory_intervals = interval_gen.generate_regulatory_intervals(gene_intervals)
    
    # Combine all intervals
    all_intervals = gene_intervals + regulatory_intervals
    
    # Write intervals files
    intervals_dir = f"{args.output_dir}/intervals"
    write_intervals_bed(all_intervals, f"{intervals_dir}/all_intervals.bed")
    write_intervals_bed(gene_intervals, f"{intervals_dir}/gene_intervals.bed")
    write_intervals_bed(regulatory_intervals, f"{intervals_dir}/regulatory_intervals.bed")
    
    # Write true peaks (enriched regions) for validation
    write_true_peaks_bed(all_intervals, f"{args.output_dir}/reference/true_peaks.bed")
    
    # Generate ChIP-seq-like FASTQ files
    print("\n--- Generating ChIP-seq-like FASTQ Files ---")
    fastq_gen = EnhancedFASTQGenerator(seed=args.seed)
    
    for sample_id in range(1, args.num_samples + 1):
        print(f"\nGenerating Sample {sample_id}...")
        
        # Use different random state for each sample
        sample_seed = args.seed + sample_id * 1000
        fastq_gen = EnhancedFASTQGenerator(seed=sample_seed)
        
        # Vary enrichment slightly between samples
        sample_enrichment = args.enrichment_factor * random.uniform(0.8, 1.2)
        
        reads = fastq_gen.generate_chip_seq_reads(
            genome_file, all_intervals, args.num_reads, args.read_length,
            enrichment_factor=sample_enrichment
        )
        
        output_fastq = f"{args.output_dir}/data/sample_{sample_id:02d}.fastq.gz"
        fastq_gen.write_fastq(reads, output_fastq)
    
    # Generate configuration files
    print("\n--- Generating Configuration Files ---")
    generate_nextflow_config(args.output_dir, chromosomes)
    
    # Generate run scripts
    run_script_content = f"""#!/bin/bash
# Enhanced Interval-based Peak Calling Pipeline Runner
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

set -e

echo "=== Enhanced Interval-based Peak Calling Pipeline ==="
echo "Directory: {args.output_dir}"
echo "Samples: {args.num_samples}"
echo "Started: $(date)"
echo

# Activate conda environment
source ~/.bashrc
conda activate biotools 2>/dev/null || echo "Note: biotools conda environment not found"

# Navigate to project directory
cd {args.output_dir}

# Build Rust project
echo "Building Rust peak caller..."
cargo build --release

# Run Nextflow pipeline
echo "Running Nextflow pipeline..."
nextflow run ../main.nf \\
    -c nextflow.config \\
    --input_dir data \\
    --output_dir results \\
    --reference_genome reference/genome.fa \\
    --intervals_bed intervals/regulatory_intervals.bed \\
    --window_size 500 \\
    --min_coverage 3.0 \\
    --pvalue_threshold 0.05 \\
    --threads 4 \\
    -with-report results/nextflow_execution_report.html \\
    -with-timeline results/nextflow_timeline.html \\
    -with-trace results/nextflow_trace.txt

echo
echo "=== Pipeline Complete ==="
echo "Results available in: {args.output_dir}/results/"
echo "Reports available in: {args.output_dir}/results/"
echo "Finished: $(date)"
"""
    
    run_script_path = f"{args.output_dir}/run_pipeline.sh"
    with open(run_script_path, 'w') as f:
        f.write(run_script_content)
    os.chmod(run_script_path, 0o755)
    
    # Generate README
    readme_content = f"""# Enhanced Interval-based Peak Calling Test Dataset

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Dataset Overview

This synthetic dataset was created for testing the enhanced interval-based peak calling pipeline.

### Genome Information
- **Chromosomes**: {len(chromosomes)}
- **Total size**: {sum(length for _, length in chromosomes):,} bp
- **Features**: Realistic GC content, repeat elements, CpG islands

### Intervals Generated
- **Gene-like intervals**: {len(gene_intervals)}
- **Regulatory intervals**: {len(regulatory_intervals)}
- **Total intervals**: {len(all_intervals)}
- **Enriched regions**: {len([i for i in all_intervals if i[4] in ['promoter', 'enhancer']])}

### Sample Information
- **Number of samples**: {args.num_samples}
- **Reads per sample**: {args.num_reads:,}
- **Read length**: {args.read_length} bp
- **Enrichment factor**: {args.enrichment_factor}x

## Directory Structure

```
{args.output_dir}/
├── data/                    # FASTQ files
├── reference/               # Reference genome and annotations
├── intervals/               # Interval definitions
├── config/                  # Configuration files
├── results/                 # Pipeline outputs (after running)
├── nextflow.config         # Pipeline configuration
├── run_pipeline.sh         # Main execution script
└── README.md               # This file
```

## Quick Start

1. **Build the pipeline**:
   ```bash
   cd {args.output_dir}
   cargo build --release
   ```

2. **Run the complete pipeline**:
   ```bash
   ./run_pipeline.sh
   ```

3. **Or run individual steps**:
   ```bash
   # Test Rust peak caller directly
   cargo run --release --bin rust_peak_caller -- \\
     --input data/sample_01_aligned.bam \\
     --intervals intervals/regulatory_intervals.bed \\
     --output test_peaks.json
   
   # Run Nextflow pipeline
   nextflow run ../main.nf -c nextflow.config
   ```

## Key Files

- `reference/genome.fa` - Synthetic reference genome
- `intervals/regulatory_intervals.bed` - Target intervals for analysis
- `reference/true_peaks.bed` - Known enriched regions for validation
- `data/sample_*.fastq.gz` - ChIP-seq-like sequencing data

## Expected Results

The pipeline should identify peaks primarily in:
- Promoter regions (upstream of genes)
- Enhancer regions (regulatory elements)

Background regions should show minimal enrichment.

## Validation

Compare called peaks with `reference/true_peaks.bed` to assess:
- Sensitivity