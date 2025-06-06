#!/usr/bin/env python3
"""
Multi-Omics Test Dataset Generator
Generates realistic multi-omics test data for integration pipeline validation
"""

import os
import gzip
import random
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime

class MultiOmicsDataGenerator:
    """Generate realistic multi-omics datasets with known correlations"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
        self.bases = ['A', 'T', 'G', 'C']
    
    def generate_reference_genome(self, chromosomes: List[Tuple[str, int]], output_file: str):
        """Generate reference genome"""
        print(f"Generating reference genome: {output_file}")
        
        with open(output_file, 'w') as f:
            for chrom_name, length in chromosomes:
                print(f"  Generating {chrom_name} ({length:,} bp)")
                f.write(f">{chrom_name}\n")
                
                # Generate sequence with realistic GC content
                sequence = self.generate_sequence(length)
                
                # Write in 80-character lines
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
    
    def generate_sequence(self, length: int) -> str:
        """Generate genomic sequence with realistic composition"""
        sequence = []
        gc_content = 0.42
        
        for i in range(length):
            # Create GC islands occasionally
            if i % 20000 < 2000 and random.random() < 0.1:
                local_gc = min(0.7, gc_content + 0.25)
            else:
                local_gc = gc_content + random.uniform(-0.1, 0.1)
            
            if random.random() < local_gc:
                sequence.append(random.choice(['G', 'C']))
            else:
                sequence.append(random.choice(['A', 'T']))
        
        return ''.join(sequence)
    
    def generate_gene_annotations(self, chromosomes: List[Tuple[str, int]], 
                                genes_per_mb: float = 25.0) -> List[Dict]:
        """Generate realistic gene annotations"""
        genes = []
        gene_types = ['protein_coding', 'lncRNA', 'miRNA', 'pseudogene', 'enhancer']
        type_weights = [0.7, 0.15, 0.05, 0.07, 0.03]
        
        for chrom_name, chrom_length in chromosomes:
            num_genes = int((chrom_length / 1_000_000) * genes_per_mb)
            
            for i in range(num_genes):
                gene_length = random.randint(1000, 100000)
                gene_start = random.randint(1000, chrom_length - gene_length - 1000)
                gene_end = gene_start + gene_length
                
                gene_type = random.choices(gene_types, weights=type_weights)[0]
                gene_id = f"GENE_{chrom_name}_{i+1:04d}"
                gene_name = f"{gene_type.upper()}_{chrom_name}_{i+1}"
                strand = random.choice(['+', '-'])
                
                genes.append({
                    'chromosome': chrom_name,
                    'start': gene_start,
                    'end': gene_end,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'strand': strand
                })
        
        return genes
    
    def write_gene_annotation_gtf(self, genes: List[Dict], output_file: str):
        """Write gene annotations in GTF format"""
        print(f"Writing {len(genes)} gene annotations to {output_file}")
        
        with open(output_file, 'w') as f:
            f.write("# GTF file generated for multi-omics testing\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            
            for gene in genes:
                # Gene line
                f.write(f"{gene['chromosome']}\ttest\tgene\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\t")
                f.write(f'gene_id "{gene['gene_id']}"; gene_name "{gene['gene_name']}"; gene_type "{gene['gene_type']}";\n')
                
                # Transcript line (simplified)
                transcript_id = f"{gene['gene_id']}_T1"
                f.write(f"{gene['chromosome']}\ttest\ttranscript\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\t")
                f.write(f'gene_id "{gene['gene_id']}"; transcript_id "{transcript_id}"; gene_name "{gene['gene_name']}"; gene_type "{gene['gene_type']}";\n')
                
                # Exon line (simplified)
                f.write(f"{gene['chromosome']}\ttest\texon\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\t")
                f.write(f'gene_id "{gene['gene_id']}"; transcript_id "{transcript_id}"; exon_number "1";\n')

class ChIPSeqGenerator:
    """Generate ChIP-seq-like BAM files with realistic binding patterns"""
    
    def __init__(self, seed: int = 42):
        random.seed(seed)
        np.random.seed(seed)
    
    def generate_chipseq_reads(self, genes: List[Dict], num_reads: int, 
                             binding_probability: float = 0.15,
                             enrichment_factor: float = 8.0) -> List[Tuple[str, int, int, str]]:
        """Generate ChIP-seq reads with enrichment at gene promoters"""
        reads = []
        
        # Select subset of genes to be "bound" by transcription factor
        bound_genes = random.sample(genes, int(len(genes) * binding_probability))
        
        enriched_reads = int(num_reads * 0.4)  # 40% from enriched regions
        background_reads = num_reads - enriched_reads
        
        print(f"Generating {enriched_reads:,} enriched ChIP-seq reads at {len(bound_genes)} bound genes")
        
        # Generate enriched reads at promoter regions
        for _ in range(enriched_reads):
            gene = random.choice(bound_genes)
            
            # Promoter region (upstream of TSS)
            if gene['strand'] == '+':
                promoter_center = gene['start'] - 1000  # 1kb upstream
            else:
                promoter_center = gene['end'] + 1000    # 1kb downstream (upstream for minus strand)
            
            # Fragment around promoter
            fragment_size = random.randint(150, 400)
            fragment_start = promoter_center - fragment_size // 2 + random.randint(-500, 500)
            fragment_end = fragment_start + fragment_size
            
            # Ensure within chromosome bounds
            fragment_start = max(1, fragment_start)
            
            read_id = f"chip_enriched_{len(reads)}"
            reads.append((gene['chromosome'], fragment_start, fragment_end, read_id))
        
        # Generate background reads
        for _ in range(background_reads):
            gene = random.choice(genes)  # Any gene for background
            
            # Random position within gene body or nearby
            gene_center = (gene['start'] + gene['end']) // 2
            fragment_size = random.randint(150, 400)
            fragment_start = gene_center + random.randint(-10000, 10000)
            fragment_end = fragment_start + fragment_size
            
            fragment_start = max(1, fragment_start)
            
            read_id = f"chip_background_{len(reads) - enriched_reads}"
            reads.append((gene['chromosome'], fragment_start, fragment_end, read_id))
        
        return reads, bound_genes
    
    def write_sam_file(self, reads: List[Tuple[str, int, int, str]], 
                      chromosomes: List[Tuple[str, int]], output_file: str):
        """Write reads to SAM format"""
        print(f"Writing {len(reads):,} ChIP-seq reads to {output_file}")
        
        with open(output_file, 'w') as f:
            # Write SAM header
            f.write("@HD\tVN:1.6\tSO:coordinate\n")
            for chrom_name, chrom_length in chromosomes:
                f.write(f"@SQ\tSN:{chrom_name}\tLN:{chrom_length}\n")
            
            # Write reads
            for chrom, start, end, read_id in reads:
                read_length = 75
                sequence = 'N' * read_length  # Placeholder sequence
                quality = 'I' * read_length   # High quality scores
                
                # SAM format: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
                f.write(f"{read_id}\t0\t{chrom}\t{start}\t60\t{read_length}M\t*\t0\t0\t{sequence}\t{quality}\n")

class ATACSeqGenerator:
    """Generate ATAC-seq