#!/usr/bin/env python3
"""
Synthetic Epigenomic Dataset Generator
Generates realistic test data for the peak calling pipeline
"""

import os
import gzip
import random
import subprocess
import argparse
from pathlib import Path
from typing import List, Tuple
import numpy as np

class GenomeGenerator:
    """Generate synthetic reference genome"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
    
    def generate_chromosome(self, length: int, gc_content: float = 0.4) -> str:
        """Generate a single chromosome sequence"""
        sequence = []
        for _ in range(length):
            if random.random() < gc_content:
                sequence.append(random.choice(['G', 'C']))
            else:
                sequence.append(random.choice(['A', 'T']))
        return ''.join(sequence)
    
    def generate_genome(self, chromosomes: List[Tuple[str, int]], output_file: str):
        """Generate multi-chromosome reference genome"""
        print(f"Generating reference genome with {len(chromosomes)} chromosomes...")
        
        with open(output_file, 'w') as f:
            for chrom_name, length in chromosomes:
                print(f"  Generating {chrom_name} ({length:,} bp)")
                f.write(f">{chrom_name}\n")
                
                # Write sequence in 80-character lines
                sequence = self.generate_chromosome(length)
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
        
        print(f"Reference genome saved to {output_file}")

class PeakRegionGenerator:
    """Generate enriched regions for peak simulation"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        self.peak_regions = []
    
    def generate_peak_regions(self, chromosomes: List[Tuple[str, int]], 
                            num_peaks_per_mb: float = 10) -> List[Tuple[str, int, int, float]]:
        """Generate random peak regions with enrichment factors"""
        peak_regions = []
        
        for chrom_name, chrom_length in chromosomes:
            # Calculate number of peaks for this chromosome
            num_peaks = int((chrom_length / 1_000_000) * num_peaks_per_mb)
            
            for _ in range(num_peaks):
                # Random peak location
                peak_center = random.randint(1000, chrom_length - 1000)
                peak_width = random.randint(200, 2000)  # Variable peak width
                peak_start = max(0, peak_center - peak_width // 2)
                peak_end = min(chrom_length, peak_center + peak_width // 2)
                
                # Random enrichment factor (2x to 20x)
                enrichment = random.uniform(2.0, 20.0)
                
                peak_regions.append((chrom_name, peak_start, peak_end, enrichment))
        
        self.peak_regions = peak_regions
        print(f"Generated {len(peak_regions)} peak regions")
        return peak_regions

class FASTQGenerator:
    """Generate synthetic FASTQ reads with realistic quality scores"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
        
    def generate_quality_string(self, length: int, base_quality: int = 30) -> str:
        """Generate realistic quality scores"""
        qualities = []
        for i in range(length):
            # Quality decreases towards the end of reads
            pos_factor = 1.0 - (i / length) * 0.3
            quality = max(20, int(base_quality * pos_factor + random.gauss(0, 3)))
            quality = min(40, quality)  # Cap at Q40
            qualities.append(chr(quality + 33))  # Convert to ASCII
        return ''.join(qualities)
    
    def mutate_sequence(self, sequence: str, error_rate: float = 0.01) -> str:
        """Introduce sequencing errors"""
        mutated = list(sequence)
        for i in range(len(mutated)):
            if random.random() < error_rate:
                mutated[i] = random.choice(self.bases)
        return ''.join(mutated)
    
    def generate_reads_from_genome(self, genome_file: str, peak_regions: List[Tuple[str, int, int, float]],
                                 num_reads: int, read_length: int = 75, 
                                 background_coverage: float = 1.0) -> List[Tuple[str, str, str]]:
        """Generate reads with enrichment in peak regions"""
        print(f"Generating {num_reads:,} reads...")
        
        # Load genome sequences
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
        
        print(f"Loaded {len(genome_seqs)} chromosomes")
        
        reads = []
        peak_reads = 0
        background_reads = 0
        
        # Create enrichment map
        enrichment_map = {}
        for chrom, start, end, enrichment in peak_regions:
            if chrom not in enrichment_map:
                enrichment_map[chrom] = []
            enrichment_map[chrom].append((start, end, enrichment))
        
        for read_id in range(num_reads):
            # Decide if this read comes from a peak region
            if random.random() < 0.3:  # 30% of reads from peaks
                # Sample from peak regions
                chrom, start, end, enrichment = random.choice(peak_regions)
                if chrom in genome_seqs:
                    # Sample position within peak region
                    max_start = max(0, end - read_length)
                    if max_start >= start:
                        read_start = random.randint(start, max_start)
                        read_seq = genome_seqs[chrom][read_start:read_start + read_length]
                        if len(read_seq) == read_length:
                            read_seq = self.mutate_sequence(read_seq)
                            quality = self.generate_quality_string(read_length, base_quality=32)
                            reads.append((f"read_{read_id}_peak", read_seq, quality))
                            peak_reads += 1
                            continue
            
            # Background read
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            if len(chrom_seq) > read_length:
                read_start = random.randint(0, len(chrom_seq) - read_length)
                read_seq = chrom_seq[read_start:read_start + read_length]
                read_seq = self.mutate_sequence(read_seq)
                quality = self.generate_quality_string(read_length, base_quality=28)
                reads.append((f"read_{read_id}_bg", read_seq, quality))
                background_reads += 1
        
        print(f"Generated {peak_reads:,} peak reads and {background_reads:,} background reads")
        return reads
    
    def write_fastq(self, reads: List[Tuple[str, str, str]], output_file: str):
        """Write reads to compressed FASTQ file"""
        print(f"Writing {len(reads):,} reads to {output_file}")
        
        with gzip.open(output_file, 'wt') as f:
            for read_id, sequence, quality in reads:
                f.write(f"@{read_id}\n")
                f.write(f"{sequence}\n")
                f.write(f"+\n")
                f.write(f"{quality}\n")

def create_directory_structure(base_dir: str):
    """Create necessary directory structure"""
    directories = [
        f"{base_dir}/data",
        f"{base_dir}/reference",
        f"{base_dir}/results"
    ]
    
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
    
    print(f"Created directory structure in {base_dir}")

def generate_adapter_file(output_file: str):
    """Generate a simple adapter file for trimming"""
    adapters = [
        ">TruSeq_Adapter_1",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        ">TruSeq_Adapter_2", 
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    ]
    
    with open(output_file, 'w') as f:
        for line in adapters:
            f.write(line + "\n")
    
    print(f"Generated adapter file: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic epigenomic dataset")
    parser.add_argument("--output-dir", default="./epigenomic_test", 
                       help="Output directory for generated data")
    parser.add_argument("--num-reads", type=int, default=500000,
                       help="Number of reads to generate")
    parser.add_argument("--num-samples", type=int, default=3,
                       help="Number of sample FASTQ files to generate")
    parser.add_argument("--read-length", type=int, default=75,
                       help="Length of generated reads")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducible data")
    
    args = parser.parse_args()
    
    print("=== Synthetic Epigenomic Dataset Generator ===")
    print(f"Output directory: {args.output_dir}")
    print(f"Number of reads per sample: {args.num_reads:,}")
    print(f"Number of samples: {args.num_samples}")
    print(f"Read length: {args.read_length}")
    print(f"Random seed: {args.seed}")
    print()
    
    # Create directory structure
    create_directory_structure(args.output_dir)
    
    # Define small chromosomes for testing
    chromosomes = [
        ("chr1", 2_000_000),  # 2 Mb
        ("chr2", 1_500_000),  # 1.5 Mb
        ("chr3", 1_000_000),  # 1 Mb
        ("chrX", 800_000),    # 0.8 Mb
    ]
    
    # Generate reference genome
    genome_gen = GenomeGenerator(seed=args.seed)
    genome_file = f"{args.output_dir}/reference/genome.fa"
    genome_gen.generate_genome(chromosomes, genome_file)
    
    # Generate peak regions
    peak_gen = PeakRegionGenerator(seed=args.seed)
    peak_regions = peak_gen.generate_peak_regions(chromosomes, num_peaks_per_mb=15)
    
    # Save peak regions for reference
    peak_file = f"{args.output_dir}/reference/true_peaks.bed"
    with open(peak_file, 'w') as f:
        f.write("# chrom\tstart\tend\tenrichment\n")
        for chrom, start, end, enrichment in peak_regions:
            f.write(f"{chrom}\t{start}\t{end}\t{enrichment:.2f}\n")
    print(f"Saved true peak regions to {peak_file}")
    
    # Generate adapter file
    adapter_file = f"{args.output_dir}/reference/adapters.fa"
    generate_adapter_file(adapter_file)
    
    # Generate FASTQ files for multiple samples
    fastq_gen = FASTQGenerator(seed=args.seed)
    
    for sample_id in range(1, args.num_samples + 1):
        print(f"\n--- Generating Sample {sample_id} ---")
        
        # Use different random state for each sample
        sample_seed = args.seed + sample_id * 1000
        fastq_gen = FASTQGenerator(seed=sample_seed)
        
        reads = fastq_gen.generate_reads_from_genome(
            genome_file, peak_regions, args.num_reads, args.read_length
        )
        
        output_fastq = f"{args.output_dir}/data/sample_{sample_id:02d}.fastq.gz"
        fastq_gen.write_fastq(reads, output_fastq)
    
    # Generate configuration files
    config_content = f"""# Epigenomic Pipeline Configuration
input_dir = "{args.output_dir}/data"
output_dir = "{args.output_dir}/results"
reference_genome = "{args.output_dir}/reference/genome.fa"
adapter_file = "{args.output_dir}/reference/adapters.fa"
window_size = 200
min_coverage = 5.0
pvalue_threshold = 0.05
threads = 4
"""
    
    with open(f"{args.output_dir}/pipeline.config", 'w') as f:
        f.write(config_content)
    
    print(f"\n=== Dataset Generation Complete ===")
    print(f"Generated files:")
    print(f"  - Reference genome: {genome_file}")
    print(f"  - True peaks: {peak_file}")
    print(f"  - Adapter sequences: {adapter_file}")
    print(f"  - {args.num_samples} FASTQ samples in {args.output_dir}/data/")
    print(f"  - Pipeline config: {args.output_dir}/pipeline.config")
    print(f"\nTo run the pipeline:")
    print(f"  cd {args.output_dir}")
    print(f"  nextflow run main.nf -c pipeline.config")

if __name__ == "__main__":
    main()