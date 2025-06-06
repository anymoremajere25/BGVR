#!/usr/bin/env python3
"""
Fragment Shift Test Dataset Generator
Generates realistic ChIP-seq test data with known fragment shifts for validation
"""

import os
import gzip
import random
import argparse
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime

class FragmentShiftDataGenerator:
    """Generate test data with specific fragment shift characteristics"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
    
    def generate_reference_genome(self, chromosomes: List[Tuple[str, int]], output_file: str):
        """Generate reference genome with realistic sequence composition"""
        print(f"Generating reference genome: {output_file}")
        
        with open(output_file, 'w') as f:
            for chrom_name, length in chromosomes:
                print(f"  Generating {chrom_name} ({length:,} bp)")
                f.write(f">{chrom_name}\n")
                
                # Generate sequence with varying GC content
                sequence = self.generate_realistic_sequence(length)
                
                # Write in 80-character lines
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
    
    def generate_realistic_sequence(self, length: int) -> str:
        """Generate genomic sequence with realistic features"""
        sequence = []
        gc_content = 0.42  # Human genome average
        
        for i in range(length):
            # Create CpG islands occasionally
            if i % 10000 < 1000 and random.random() < 0.1:  # CpG island
                local_gc = min(0.7, gc_content + 0.25)
            # Create AT-rich regions
            elif i % 8000 < 500 and random.random() < 0.1:  # AT-rich
                local_gc = max(0.25, gc_content - 0.15)
            else:
                local_gc = gc_content + random.uniform(-0.1, 0.1)
            
            if random.random() < local_gc:
                sequence.append(random.choice(['G', 'C']))
            else:
                sequence.append(random.choice(['A', 'T']))
        
        return ''.join(sequence)
    
    def generate_binding_sites(self, chromosomes: List[Tuple[str, int]], 
                             sites_per_mb: float = 10.0) -> List[Tuple[str, int, int, float]]:
        """Generate transcription factor binding sites with realistic distribution"""
        binding_sites = []
        
        for chrom_name, chrom_length in chromosomes:
            num_sites = int((chrom_length / 1_000_000) * sites_per_mb)
            
            for i in range(num_sites):
                # Binding site width (typically narrow for TF ChIP-seq)
                site_width = random.randint(50, 500)
                site_center = random.randint(site_width, chrom_length - site_width)
                site_start = site_center - site_width // 2
                site_end = site_center + site_width // 2
                
                # Binding strength (affects fragment density)
                binding_strength = random.uniform(2.0, 20.0)
                
                binding_sites.append((chrom_name, site_start, site_end, binding_strength))
        
        return binding_sites

class ChIPSeqReadGenerator:
    """Generate ChIP-seq reads with realistic fragment characteristics"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
    
    def generate_chipseq_fragments(self, genome_seqs: Dict[str, str], 
                                 binding_sites: List[Tuple[str, int, int, float]],
                                 target_fragment_size: int,
                                 fragment_size_std: int,
                                 num_reads: int,
                                 read_length: int = 75) -> List[Tuple[str, str, str, bool, int]]:
        """Generate ChIP-seq fragments with specific size distribution"""
        fragments = []
        
        # Calculate enrichment vs background ratio
        enriched_fraction = 0.3  # 30% of reads from binding sites
        enriched_reads = int(num_reads * enriched_fraction)
        background_reads = num_reads - enriched_reads
        
        print(f"Generating {enriched_reads:,} enriched fragments and {background_reads:,} background fragments")
        
        # Generate enriched fragments at binding sites
        for i in range(enriched_reads):
            if not binding_sites:
                break
                
            chrom, site_start, site_end, strength = random.choice(binding_sites)
            if chrom not in genome_seqs:
                continue
            
            # Sample fragment size from normal distribution
            fragment_size = max(read_length + 10, 
                              int(np.random.normal(target_fragment_size, fragment_size_std)))
            
            # Position fragment around binding site
            site_center = (site_start + site_end) // 2
            fragment_start = site_center - fragment_size // 2 + random.randint(-100, 100)
            fragment_start = max(0, fragment_start)
            fragment_end = fragment_start + fragment_size
            
            if fragment_end <= len(genome_seqs[chrom]):
                fragments.append(self.create_fragment_reads(
                    chrom, fragment_start, fragment_end, 
                    genome_seqs[chrom], read_length, i, True
                ))
        
        # Generate background fragments
        for i in range(background_reads):
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            
            # Background fragments - more variable size
            fragment_size = max(read_length + 10,
                              int(np.random.normal(target_fragment_size * 0.8, fragment_size_std * 1.5)))
            
            if len(chrom_seq) > fragment_size + 1000:
                fragment_start = random.randint(500, len(chrom_seq) - fragment_size - 500)
                fragment_end = fragment_start + fragment_size
                
                fragments.append(self.create_fragment_reads(
                    chrom, fragment_start, fragment_end,
                    chrom_seq, read_length, enriched_reads + i, False
                ))
        
        return [f for f in fragments if f is not None]
    
    def create_fragment_reads(self, chrom: str, start: int, end: int, 
                            chrom_seq: str, read_length: int, 
                            fragment_id: int, is_enriched: bool) -> Tuple[str, str, str, bool, int]:
        """Create paired-end reads from a fragment"""
        fragment_size = end - start
        
        # Extract fragment sequence
        if end > len(chrom_seq):
            return None
            
        fragment_seq = chrom_seq[start:end]
        
        # Create forward read (5' end)
        forward_seq = fragment_seq[:read_length]
        forward_seq = self.add_sequencing_errors(forward_seq, 0.005 if is_enriched else 0.01)
        forward_qual = self.generate_quality_scores(len(forward_seq), 32 if is_enriched else 28)
        
        # Create reverse read (3' end, reverse complement)
        reverse_seq = fragment_seq[-read_length:]
        reverse_seq = self.reverse_complement(reverse_seq)
        reverse_seq = self.add_sequencing_errors(reverse_seq, 0.005 if is_enriched else 0.01)
        reverse_qual = self.generate_quality_scores(len(reverse_seq), 32 if is_enriched else 28)
        
        return [
            (f"fragment_{fragment_id}_R1", forward_seq, forward_qual, False, start),
            (f"fragment_{fragment_id}_R2", reverse_seq, reverse_qual, True, end - read_length)
        ]
    
    def reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def add_sequencing_errors(self, sequence: str, error_rate: float = 0.01) -> str:
        """Add realistic sequencing errors"""
        mutated = list(sequence)
        for i in range(len(mutated)):
            if random.random() < error_rate:
                mutated[i] = random.choice([b for b in self.bases if b != mutated[i]])
        return ''.join(mutated)
    
    def generate_quality_scores(self, length: int, base_quality: int = 30) -> str:
        """Generate realistic quality scores with position-dependent degradation"""
        qualities = []
        for i in range(length):
            # Quality decreases towards end with some noise
            pos_factor = 1.0 - (i / length) * 0.35
            cycle_effect = 0.05 * np.sin(i * 0.2)  # Sequencing cycle effects
            random_noise = random.gauss(0, 1.5)
            
            quality = int(base_quality * pos_factor + cycle_effect + random_noise)
            quality = max(15, min(40, quality))
            qualities.append(chr(quality + 33))
        
        return ''.join(qualities)
    
    def convert_to_single_end_reads(self, fragments: List, target_shift: int) -> List[Tuple[str, str, str]]:
        """Convert paired-end fragments to single-end reads simulating ChIP-seq"""
        single_reads = []
        
        for fragment_reads in fragments:
            if not fragment_reads or len(fragment_reads) != 2:
                continue
                
            forward_read, reverse_read = fragment_reads
            
            # Randomly select which end to sequence (simulating single-end ChIP-seq)
            if random.random() < 0.5:
                # Use forward read
                read_id, sequence, quality, is_reverse, pos = forward_read
                single_reads.append((read_id + "_fwd", sequence, quality))
            else:
                # Use reverse read  
                read_id, sequence, quality, is_reverse, pos = reverse_read
                single_reads.append((read_id + "_rev", sequence, quality))
        
        return single_reads
    
    def write_fastq(self, reads: List[Tuple[str, str, str]], output_file: str):
        """Write reads to compressed FASTQ file"""
        print(f"Writing {len(reads):,} reads to {output_file}")
        
        with gzip.open(output_file, 'wt') as f:
            for read_id, sequence, quality in reads:
                f.write(f"@{read_id}\n")
                f.write(f"{sequence}\n")
                f.write(f"+\n")
                f.write(f"{quality}\n")

def load_genome_sequences(genome_file: str) -> Dict[str, str]:
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
    
    print(f"Loaded {len(genome_seqs)} chromosomes from genome file")
    return genome_seqs

def write_bed_file(intervals: List[Tuple], output_file: str, header_comment: str):
    """Write intervals to BED format"""
    with open(output_file, 'w') as f:
        f.write(f"# {header_comment}\n")
        for interval in intervals:
            if len(interval) == 4:
                f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\n")
            elif len(interval) == 5:
                f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\t{interval[4]}\n")
    
    print(f"Wrote {len(intervals)} intervals to {output_file}")

def create_directory_structure(base_dir: str):
    """Create directory structure for shift testing"""
    directories = [
        f"{base_dir}/data",
        f"{base_dir}/reference",
        f"{base_dir}/results",
        f"{base_dir}/analysis",
        f"{base_dir}/config"
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
    
    print(f"Created directory structure in {base_dir}")

def generate_nextflow_config(base_dir: str, fragment_shifts: Dict[str, int]):
    """Generate Nextflow configuration for shift testing"""
    shift_info = ", ".join([f"{sample}: {shift}bp" for sample, shift in fragment_shifts.items()])
    
    config_content = f"""// Fragment Shift Test Configuration
// Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
// True fragment shifts: {shift_info}

params {{
    // Input/Output
    input_dir = "{base_dir}/data"
    output_dir = "{base_dir}/results"
    reference_genome = "{base_dir}/reference/genome.fa"
    
    // Shift estimation parameters
    max_shift = 600
    min_mapq = 10
    bin_size = 1
    use_fft = true
    smoothing_window = 5
    sampling_factor = 1
    
    // Peak calling parameters
    window_size = 200
    min_coverage = 3.0
    pvalue_threshold = 0.05
    min_peak_width = 50
    max_peak_width = 3000
    
    // Analysis scope
    target_chromosome = null  // Analyze all chromosomes
    analysis_start = null
    analysis_end = null
    
    // Performance
    threads = 4
    memory = '8 GB'
    time_limit = '4h'
}}

// Process-specific configurations
process {{
    withName: ESTIMATE_FRAGMENT_SHIFT {{
        cpus = 4
        memory = '8 GB'
        time = '2h'
    }}
    
    withName: SHIFT_AWARE_PEAK_CALLING {{
        cpus = 4
        memory = '6 GB'
        time = '3h'
    }}
}}

// Expected results for validation:
{chr(10).join([f"// {sample}: ~{shift} bp shift" for sample, shift in fragment_shifts.items()])}
"""
    
    config_file = f"{base_dir}/nextflow.config"
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    print(f"Generated Nextflow configuration: {config_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate fragment shift test dataset")
    parser.add_argument("--output-dir", default="./shift_test", 
                       help="Output directory for test data")
    parser.add_argument("--num-reads", type=int, default=200000,
                       help="Number of reads per sample")
    parser.add_argument("--read-length", type=int, default=75,
                       help="Length of sequencing reads")
    parser.add_argument("--fragment-shifts", nargs='+', type=int,
                       default=[100, 150, 200],
                       help="Fragment shifts to generate (bp)")
    parser.add_argument("--fragment-size-std", type=int, default=30,
                       help="Standard deviation of fragment size")
    parser.add_argument("--binding-sites-per-mb", type=float, default=15.0,
                       help="Number of binding sites per megabase")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducible data")
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("Fragment Shift Test Dataset Generator")
    print("=" * 60)
    print(f"Output directory: {args.output_dir}")
    print(f"Reads per sample: {args.num_reads:,}")
    print(f"Read length: {args.read_length} bp")
    print(f"Fragment shifts: {args.fragment_shifts} bp")
    print(f"Fragment size std: {args.fragment_size_std} bp")
    print(f"Binding sites per Mb: {args.binding_sites_per_mb}")
    print(f"Random seed: {args.seed}")
    print()
    
    # Create directory structure
    create_directory_structure(args.output_dir)
    
    # Define test chromosomes (smaller for testing)
    chromosomes = [
        ("chr1", 3_000_000),   # 3 Mb
        ("chr2", 2_500_000),   # 2.5 Mb
        ("chr3", 2_000_000),   # 2 Mb
        ("chr4", 1_500_000),   # 1.5 Mb
        ("chrX", 1_200_000),   # 1.2 Mb
    ]
    
    print("\n--- Generating Reference Genome ---")
    data_gen = FragmentShiftDataGenerator(seed=args.seed)
    genome_file = f"{args.output_dir}/reference/genome.fa"
    data_gen.generate_reference_genome(chromosomes, genome_file)
    
    print("\n--- Generating Binding Sites ---")
    binding_sites = data_gen.generate_binding_sites(chromosomes, args.binding_sites_per_mb)
    write_bed_file(binding_sites, f"{args.output_dir}/reference/binding_sites.bed",
                   "Transcription factor binding sites with strength scores")
    
    print("\n--- Loading Genome Sequences ---")
    genome_seqs = load_genome_sequences(genome_file)
    
    print("\n--- Generating ChIP-seq Samples with Different Fragment Shifts ---")
    read_gen = ChIPSeqReadGenerator(seed=args.seed)
    fragment_shifts = {}
    
    for i, target_shift in enumerate(args.fragment_shifts):
        sample_name = f"sample_shift_{target_shift}bp"
        fragment_shifts[sample_name] = target_shift
        
        print(f"\nGenerating {sample_name} (target shift: {target_shift} bp)...")
        
        # Use different random state for each sample
        sample_seed = args.seed + i * 1000
        sample_read_gen = ChIPSeqReadGenerator(seed=sample_seed)
        
        # Generate fragments with specific size distribution
        fragments = sample_read_gen.generate_chipseq_fragments(
            genome_seqs, binding_sites, target_shift, 
            args.fragment_size_std, args.num_reads, args.read_length
        )
        
        # Convert to single-end reads
        single_reads = sample_read_gen.convert_to_single_end_reads(fragments, target_shift)
        
        # Write FASTQ file
        output_fastq = f"{args.output_dir}/data/{sample_name}.fastq.gz"
        sample_read_gen.write_fastq(single_reads, output_fastq)
        
        # Generate sample description
        sample_desc = f"""Sample: {sample_name}
Target Fragment Shift: {target_shift} bp
Fragment Size Std: {args.fragment_size_std} bp
Generated Reads: {len(single_reads):,}
Read Length: {args.read_length} bp
Binding Sites: {len(binding_sites)}
Expected Cross-correlation Peak: Around {target_shift} bp
"""
        
        with open(f"{args.output_dir}/data/{sample_name}_description.txt", 'w') as f:
            f.write(sample_desc)
    
    print("\n--- Generating Configuration Files ---")
    generate_nextflow_config(args.output_dir, fragment_shifts)
    
    # Generate validation script
    validation_script = f"""#!/bin/bash
# Fragment Shift Validation Script
# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

set -e

echo "=== Fragment Shift Estimation Validation ==="
echo "Directory: {args.output_dir}"
echo "Expected shifts: {', '.join([f'{s}bp' for s in args.fragment_shifts])}"
echo "Started: $(date)"
echo

# Activate conda environment
source ~/.bashrc
conda activate biotools 2>/dev/null || echo "Note: biotools environment not found"

# Navigate to test directory
cd {args.output_dir}

# Build Rust shift estimator
echo "Building Rust tools..."
cargo build --release

# Run Nextflow pipeline
echo "Running fragment shift estimation pipeline..."
nextflow run ../main.nf \\
    -c nextflow.config \\
    --input_dir data \\
    --output_dir results \\
    --reference_genome reference/genome.fa \\
    --use_fft true \\
    --threads 4 \\
    -with-report results/execution_report.html \\
    -with-timeline results/timeline.html \\
    -with-trace results/trace.txt

echo
echo "=== Validation Results ==="
echo "Comparing estimated vs true fragment shifts:"

# Extract estimated shifts and compare with true values
for sample in data/*.fastq.gz; do
    sample_name=$(basename $sample .fastq.gz)
    true_shift=$(echo $sample_name | grep -o '[0-9]\\+')
    
    if [ -f "results/shift_analysis/${{sample_name}}_shift_estimate.json" ]; then
        estimated_shift=$(jq -r '.estimated_shift' "results/shift_analysis/${{sample_name}}_shift_estimate.json")
        confidence=$(jq -r '.confidence_score' "results/shift_analysis/${{sample_name}}_shift_estimate.json")
        
        echo "$sample_name:"
        echo "  True shift: ${{true_shift}} bp"
        echo "  Estimated shift: ${{estimated_shift}} bp"
        echo "  Confidence: $confidence"
        echo "  Accuracy: $(echo "scale=1; 100 - 100 * sqrt(($estimated_shift - $true_shift)^2) / $true_shift" | bc -l)%"
        echo
    else
        echo "$sample_name: ERROR - No shift estimate found"
    fi
done

echo "Finished: $(date)"
echo "Results available in: {args.output_dir}/results/"
echo "HTML report: {args.output_dir}/results/summary/peak_calling_report.html"
"""
    
    validation_script_path = f"{args.output_dir}/run_validation.sh"
    with open(validation_script_path, 'w') as f:
        f.write(validation_script)
    os.chmod(validation_script_path, 0o755)
    
    # Generate README
    readme_content = f"""# Fragment Shift Estimation Test Dataset

Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Overview

This test dataset contains ChIP-seq samples with known fragment shifts to validate the cross-correlation-based shift estimation pipeline.

## Dataset Information

- **Genome size**: {sum(length for _, length in chromosomes):,} bp ({len(chromosomes)} chromosomes)
- **Binding sites**: {len(binding_sites)} sites
- **Reads per sample**: {args.num_reads:,}
- **Read length**: {args.read_length} bp
- **Fragment shifts tested**: {args.fragment_shifts} bp

## Fragment Shift Scenarios

"""
    
    for shift in args.fragment_shifts:
        readme_content += f"""
### sample_shift_{shift}bp.fastq.gz
- **True fragment shift**: {shift} bp
- **Expected cross-correlation peak**: Around {shift} bp
- **Fragment size distribution**: {shift} ± {args.fragment_size_std} bp
- **Validation**: Estimated shift should be within 10-20 bp of true value
"""
    
    readme_content += f"""

## Directory Structure

```
{args.output_dir}/
├── data/                    # FASTQ files with different fragment shifts
├── reference/               # Reference genome and binding sites
├── results/                 # Pipeline outputs (after running)
├── nextflow.config         # Pipeline configuration
├── run_validation.sh       # Validation script
└── README.md               # This file
```

## Quick Start

1. **Build tools**:
   ```bash
   cd {args.output_dir}
   cargo build --release
   ```

2. **Run validation**:
   ```bash
   ./run_validation.sh
   ```

3. **Manual analysis** (single sample):
   ```bash
   # Test shift estimator directly
   cargo run --release --bin shift_estimator -- \\
     --input results/alignments/sample_shift_150bp_aligned.bam \\
     --output test_shift.json \\
     --max-shift 400 \\
     --use-fft \\
     --verbose
   
   # Check estimated shift
   jq '.estimated_shift' test_shift.json
   ```

## Expected Results

The pipeline should accurately estimate fragment shifts:

| Sample | True Shift | Expected Accuracy |
|--------|------------|-------------------|
{chr(10).join([f"| sample_shift_{s}bp | {s} bp | ±10-20 bp |" for s in args.fragment_shifts])}

## Key Validation Metrics

1. **Shift Accuracy**: Estimated shift within 15% of true value
2. **Confidence Score**: Should be > 0.6 for good quality data
3. **Signal-to-Noise Ratio**: Should be > 3.0 for clear peaks
4. **Cross-correlation Profile**: Clear peak at expected shift

## Output Files

- `results/shift_analysis/` - Individual shift estimates
- `results/analysis/shift_comparison.txt` - Comparison across samples
- `results/summary/peak_calling_report.html` - Comprehensive report
- `results/peaks/` - Shift-corrected peak calls

## Troubleshooting

- **Low confidence scores**: Check data quality and increase read count
- **Inaccurate shifts**: Verify fragment size distribution and binding site density
- **Missing peaks**: Adjust cross-correlation parameters or max shift range

## Advanced Usage

Test different parameters:
```bash
# Test without FFT acceleration
nextflow run ../main.nf -c nextflow.config --use_fft false

# Test with different max shift range
nextflow run ../main.nf -c nextflow.config --max_shift 800

# Test chromosome-specific analysis
nextflow run ../main.nf -c nextflow.config --target_chromosome chr1
```

This dataset provides ground truth for validating fragment shift estimation algorithms and optimizing cross-correlation parameters.
"""
    
    readme_path = f"{args.output_dir}/README.md"
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    
    print(f"\n{'=' * 60}")
    print("🎉 Fragment Shift Test Dataset Generation Complete!")
    print(f"{'=' * 60}")
    print(f"📁 Dataset location: {args.output_dir}")
    print(f"🧬 Reference genome: {genome_file}")
    print(f"🎯 Binding sites: {len(binding_sites)} targets")
    print(f"📊 Fragment shifts: {len(args.fragment_shifts)} scenarios")
    print(f"📈 Expected accuracy: ±10-20 bp")
    print()
    print("📝 Next steps:")
    print(f"   1. cd {args.output_dir}")
    print("   2. Copy main.rs, peak_caller.rs, main.nf, and Cargo.toml")
    print("   3. ./run_validation.sh")
    print()
    print("🔬 Expected validation:")
    for shift in args.fragment_shifts:
        print(f"   • sample_shift_{shift}bp: ~{shift} bp estimated")

if __name__ == "__main__":
    main()