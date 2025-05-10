#!/usr/bin/env nextflow

// Define parameters with defaults
params.reference = 'reference.fasta'
params.region_list = 'regions.txt'
params.output_dir = 'results'
params.threads = Runtime.runtime.availableProcessors()
params.memory = '2.GB'
params.container_version = 'latest'

log.info """
🧬 GENOMIC ANALYSIS PIPELINE
Reference: ${params.reference}
Regions: ${params.region_list}
Output: ${params.output_dir}
"""

// Create input channels
reference_file = file(params.reference)
regionChannel = Channel.fromPath(params.region_list)
    .splitText()
    .map { it.trim() }
    .filter { it.length() > 0 }

// Create a simple tool process
process createTool {
    output:
    path "dummy_tool.sh", emit: tool

    script:
    """
    cat > dummy_tool.sh << 'EOF'
#!/bin/bash
# Simple dummy tool that simulates the behavior of rust_mmap_tool
reference=\$2
region=\$4
threads=\$6
output=\$8

echo "Processing region \$region with \$threads threads" >&2
echo "Reference: \$reference" >&2

# Generate dummy coverage data
echo "\$region 10" > \$output
echo "\$region 15" >> \$output
echo "\$region 20" >> \$output

exit 0
EOF
    chmod +x dummy_tool.sh
    """
}

// Coverage mapping process
process coverageMapping {
    publishDir "${params.output_dir}/coverage", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 3
    cpus params.threads
    memory params.memory

    input:
    val region
    path reference
    path tool

    output:
    path "coverage_${region.replaceAll(':', '_')}.txt"

    script:
    """
    ./${tool} \\
        --reference ${reference} \\
        --region ${region} \\
        --threads ${task.cpus} \\
        --output coverage_${region.replaceAll(':', '_')}.txt
    """
}

// Merge coverage results
process mergeCoverage {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path coverage_files
    
    output:
    path "merged_coverage.txt"
    path "coverage_summary.json"
    
    script:
    """
    cat ${coverage_files.join(" ")} > merged_coverage.txt
    
    # Simple Python script to generate summary
    python3 -c "
import json
import statistics

counts = []
with open('merged_coverage.txt', 'r') as f:
    for line in f:
        if line.strip():
            parts = line.strip().split()
            if len(parts) > 1:
                try:
                    count = int(parts[-1])
                    counts.append(count)
                except ValueError:
                    pass

summary = {
    'total_regions': len(counts),
    'min_coverage': min(counts) if counts else 0,
    'max_coverage': max(counts) if counts else 0,
    'mean_coverage': statistics.mean(counts) if counts else 0,
    'median_coverage': statistics.median(counts) if counts else 0
}

with open('coverage_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
    " || echo '{"error": "Failed to generate summary"}' > coverage_summary.json
    """
}

// Define the workflow
workflow {
    createTool()
    coverage_results = coverageMapping(regionChannel, reference_file, createTool.out.tool)
    mergeCoverage(coverage_results.collect())
}

// Completion handler
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}, duration: ${workflow.duration}"
    if (!workflow.success) {
        log.error "Pipeline execution failed"
    }
}