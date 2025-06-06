#!/usr/bin/env python3
"""
Comprehensive QC Test Dataset Generator
Generates realistic test data with various quality scenarios for testing the QC pipeline
"""

import os
import gzip
import random
import argparse
import numpy as np
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime

class QCTestDataGenerator:
    """Generate test datasets with different quality characteristics for QC validation"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
        
    def generate_reference_genome(self, chromosomes: List[Tuple[str, int]], output_file: str):
        """Generate reference genome with realistic features"""
        print(f"Generating reference genome: {output_file}")
        
        with open(output_file, 'w') as f:
            for chrom_name, length in chromosomes:
                print(f"  Generating {chrom_name} ({length:,} bp)")
                f.write(f">{chrom_name}\n")
                
                # Generate sequence with varying GC content
                gc_content = 0.41 + random.uniform(-0.05, 0.05)
                sequence = self.generate_sequence_with_features(length, gc_content)
                
                # Write in 80-character lines
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
    
    def generate_sequence_with_features(self, length: int, gc_content: float) -> str:
        """Generate sequence with realistic genomic features"""
        sequence = []
        
        for i in range(length):
            # Create occasional GC-rich regions (CpG islands)
            if i % 5000 < 500 and random.random() < 0.2:  # GC island
                local_gc = min(0.8, gc_content + 0.3)
            # Create AT-rich regions
            elif i % 3000 < 200 and random.random() < 0.1:  # AT-rich region
                local_gc = max(0.2, gc_content - 0.2)
            else:
                local_gc = gc_content
                
            if random.random() < local_gc:
                sequence.append(random.choice(['G', 'C']))
            else:
                sequence.append(random.choice(['A', 'T']))
        
        return ''.join(sequence)
    
    def generate_peak_regions(self, chromosomes: List[Tuple[str, int]], 
                            num_peaks_per_mb: float = 15) -> List[Tuple[str, int, int, str, float]]:
        """Generate realistic peak regions with varying enrichment"""
        peak_regions = []
        
        for chrom_name, chrom_length in chromosomes:
            num_peaks = int((chrom_length / 1_000_000) * num_peaks_per_mb)
            
            for i in range(num_peaks):
                # Varying peak sizes
                peak_width = random.randint(300, 3000)
                peak_center = random.randint(peak_width, chrom_length - peak_width)
                peak_start = peak_center - peak_width // 2
                peak_end = peak_center + peak_width // 2
                
                # Varying enrichment factors
                enrichment = random.uniform(3.0, 25.0)
                peak_name = f"peak_{chrom_name}_{i+1}"
                
                peak_regions.append((chrom_name, peak_start, peak_end, peak_name, enrichment))
        
        return peak_regions
    
    def generate_background_regions(self, chromosomes: List[Tuple[str, int]], 
                                  num_regions: int = 500) -> List[Tuple[str, int, int, str]]:
        """Generate background regions for QC analysis"""
        background_regions = []
        
        for _ in range(num_regions):
            chrom_name, chrom_length = random.choice(chromosomes)
            region_size = 1000  # 1kb background regions
            
            if chrom_length > region_size * 2:
                start = random.randint(1000, chrom_length - region_size - 1000)
                end = start + region_size
                name = f"background_{len(background_regions) + 1}"
                
                background_regions.append((chrom_name, start, end, name))
        
        return background_regions
    
    def generate_blacklist_regions(self, chromosomes: List[Tuple[str, int]], 
                                 num_regions: int = 50) -> List[Tuple[str, int, int, str]]:
        """Generate blacklist regions (repetitive/problematic regions)"""
        blacklist_regions = []
        
        for _ in range(num_regions):
            chrom_name, chrom_length = random.choice(chromosomes)
            region_size = random.randint(500, 5000)  # Variable size blacklist regions
            
            if chrom_length > region_size * 2:
                start = random.randint(0, chrom_length - region_size)
                end = start + region_size
                name = f"blacklist_{len(blacklist_regions) + 1}"
                
                blacklist_regions.append((chrom_name, start, end, name))
        
        return blacklist_regions

class QCScenarioGenerator:
    """Generate FASTQ files with different quality scenarios"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
    
    def generate_high_quality_sample(self, genome_seqs: Dict[str, str], peak_regions: List,
                                   num_reads: int, read_length: int = 75) -> List[Tuple[str, str, str]]:
        """Generate high-quality ChIP-seq sample (good enrichment, low duplication)"""
        reads = []
        
        # 35% reads from peaks, 65% background
        peak_reads = int(num_reads * 0.35)
        background_reads = num_reads - peak_reads
        
        # Generate peak reads with good enrichment
        for i in range(peak_reads):
            chrom, start, end, name, enrichment = random.choice(peak_regions)
            if chrom in genome_seqs:
                read_start = random.randint(start, max(start, end - read_length))
                read_seq = genome_seqs[chrom][read_start:read_start + read_length]
                read_seq = self.add_sequencing_errors(read_seq, error_rate=0.005)  # Low error rate
                quality = self.generate_quality_scores(read_length, base_quality=35)  # High quality
                reads.append((f"read_{i}_peak_hq", read_seq, quality))
        
        # Generate background reads
        for i in range(background_reads):
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            read_start = random.randint(0, len(chrom_seq) - read_length)
            read_seq = chrom_seq[read_start:read_start + read_length]
            read_seq = self.add_sequencing_errors(read_seq, error_rate=0.008)
            quality = self.generate_quality_scores(read_length, base_quality=32)
            reads.append((f"read_{i}_bg_hq", read_seq, quality))
        
        return reads
    
    def generate_low_quality_sample(self, genome_seqs: Dict[str, str], peak_regions: List,
                                  num_reads: int, read_length: int = 75) -> List[Tuple[str, str, str]]:
        """Generate low-quality sample (poor enrichment, high duplication, low mapping)"""
        reads = []
        
        # Only 5% reads from peaks (poor enrichment)
        peak_reads = int(num_reads * 0.05)
        background_reads = num_reads - peak_reads
        
        # Generate few peak reads
        for i in range(peak_reads):
            chrom, start, end, name, enrichment = random.choice(peak_regions)
            if chrom in genome_seqs:
                read_start = random.randint(start, max(start, end - read_length))
                read_seq = genome_seqs[chrom][read_start:read_start + read_length]
                read_seq = self.add_sequencing_errors(read_seq, error_rate=0.02)  # High error rate
                quality = self.generate_quality_scores(read_length, base_quality=25)  # Lower quality
                reads.append((f"read_{i}_peak_lq", read_seq, quality))
        
        # Generate background reads with high duplication
        unique_background = background_reads // 3  # High duplication rate
        for i in range(unique_background):
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            read_start = random.randint(0, len(chrom_seq) - read_length)
            read_seq = chrom_seq[read_start:read_start + read_length]
            read_seq = self.add_sequencing_errors(read_seq, error_rate=0.015)
            quality = self.generate_quality_scores(read_length, base_quality=28)
            
            # Create duplicates
            dup_count = random.randint(2, 4)
            for d in range(dup_count):
                reads.append((f"read_{i}_bg_lq_dup{d}", read_seq, quality))
                if len(reads) >= num_reads:
                    break
            if len(reads) >= num_reads:
                break
        
        # Add some random/unmappable sequences
        random_reads = min(num_reads - len(reads), num_reads // 10)
        for i in range(random_reads):
            random_seq = ''.join(random.choices(self.bases, k=read_length))
            quality = self.generate_quality_scores(read_length, base_quality=20)
            reads.append((f"read_{i}_random", random_seq, quality))
        
        return reads[:num_reads]
    
    def generate_medium_quality_sample(self, genome_seqs: Dict[str, str], peak_regions: List,
                                     num_reads: int, read_length: int = 75) -> List[Tuple[str, str, str]]:
        """Generate medium-quality sample (moderate enrichment and duplication)"""
        reads = []
        
        # 15% reads from peaks (moderate enrichment)
        peak_reads = int(num_reads * 0.15)
        background_reads = num_reads - peak_reads
        
        # Generate peak reads
        for i in range(peak_reads):
            chrom, start, end, name, enrichment = random.choice(peak_regions)
            if chrom in genome_seqs:
                read_start = random.randint(start, max(start, end - read_length))
                read_seq = genome_seqs[chrom][read_start:read_start + read_length]
                read_seq = self.add_sequencing_errors(read_seq, error_rate=0.01)
                quality = self.generate_quality_scores(read_length, base_quality=30)
                reads.append((f"read_{i}_peak_mq", read_seq, quality))
        
        # Generate background reads with moderate duplication
        unique_background = int(background_reads * 0.7)  # 30% duplication
        for i in range(unique_background):
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            read_start = random.randint(0, len(chrom_seq) - read_length)
            read_seq = chrom_seq[read_start:read_start + read_length]
            read_seq = self.add_sequencing_errors(read_seq, error_rate=0.01)
            quality = self.generate_quality_scores(read_length, base_quality=30)
            
            reads.append((f"read_{i}_bg_mq", read_seq, quality))
            
            # Some duplicates
            if random.random() < 0.3:  # 30% chance of duplication
                reads.append((f"read_{i}_bg_mq_dup", read_seq, quality))
        
        # Fill remaining with unique reads
        while len(reads) < num_reads:
            chrom = random.choice(list(genome_seqs.keys()))
            chrom_seq = genome_seqs[chrom]
            read_start = random.randint(0, len(chrom_seq) - read_length)
            read_seq = chrom_seq[read_start:read_start + read_length]
            read_seq = self.add_sequencing_errors(read_seq, error_rate=0.01)
            quality = self.generate_quality_scores(read_length, base_quality=30)
            reads.append((f"read_{len(reads)}_bg_mq_fill", read_seq, quality))
        
        return reads[:num_reads]
    
    def add_sequencing_errors(self, sequence: str, error_rate: float = 0.01) -> str:
        """Add sequencing errors to reads"""
        mutated = list(sequence)
        for i in range(len(mutated)):
            if random.random() < error_rate:
                # Random substitution
                mutated[i] = random.choice([b for b in self.bases if b != mutated[i]])
        return ''.join(mutated)
    
    def generate_quality_scores(self, length: int, base_quality: int = 30) -> str:
        """Generate realistic quality scores"""
        qualities = []
        for i in range(length):
            # Quality decreases towards end
            pos_factor = 1.0 - (i / length) * 0.3
            noise = random.gauss(0, 2)
            quality = int(base_quality * pos_factor + noise)
            quality = max(15, min(40, quality))
            qualities.append(chr(quality + 33))
        return ''.join(qualities)
    
    def write_fastq(self, reads: List[Tuple[str, str, str]], output_file: str):
        """Write reads to compressed FASTQ file"""
        print(f"Writing {len(reads):,} reads to {output_file}")
        
        with gzip.open(output_file, 'wt') as f:
            for read_id, sequence, quality in reads:
                f.write(f"@{read_id}\n")
                f.write(f"{sequence}\n")
                f.write(f"+\n")
                f.write(f"{quality}\n")

def write_bed_file(intervals: List[Tuple], output_file: str, header_comment: str):
    """Write intervals to BED format"""
    with open(output_file, 'w') as f:
        f.write(f"# {header_comment}\n")
        for interval in intervals:
            if len(interval) == 4:  # (chrom, start, end, name)
                f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\n")
            elif len(interval) == 5:  # (chrom, start, end, name, score)
                f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\t{interval[3]}\t{interval[4]}\n")
    
    print(f"Wrote {len(intervals)} intervals to {output_file}")

def load_genome_sequences(genome_file: str) -> Dict[str, str]:
    """Load genome sequences from FASTA file"""
    genome