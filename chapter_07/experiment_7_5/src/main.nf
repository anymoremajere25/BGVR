#!/usr/bin/env nextflow

//
// A minimal Nextflow script to run a local Rust coverage tool,
// parallelizing tasks over multiple BAM files and merging the JSON outputs.
//

process RUN_COVERAGE {
    // Use the local system instead of a container
    container = null

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)

    output:
    path("${sample_id}.coverage.json")

    script:
    """
    echo "Running coverage on sample: ${sample_id}"
    /mnt/c/Users/trian/BGVR/chapter_07/experiment_75/target/debug/rust_coverage_tool --bam ${bam_file} --bai ${bai_file} --region chr1:10000-10100 > ${sample_id}.coverage.json
    """
}

process MERGE_COVERAGE {
    // Use the local system without jq dependency
    container = null

    input:
    path coverage_files

    output:
    path('merged_coverage.json')

    script:
    """
    # Simple merge of JSON files using bash
    echo "[" > merged_coverage.json
    
    # Add each file content with commas between
    COUNT=\$(ls -1 ${coverage_files} | wc -l)
    CURRENT=0
    
    for file in ${coverage_files}; do
        CURRENT=\$((CURRENT+1))
        cat \$file >> merged_coverage.json
        
        # Add comma if not the last file
        if [ \$CURRENT -lt \$COUNT ]; then
            echo "," >> merged_coverage.json
        fi
    done
    
    echo "]" >> merged_coverage.json
    """
}

workflow {
    // Channel with explicit BAM and BAI file pairs
    bam_files = Channel.fromPath('*.bam')
        .map { bam ->
            def sample_id = bam.simpleName
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                bai = file("${bam.parent}/${sample_id}.bai")
            }
            tuple(sample_id, bam, bai)
        }

    // Run coverage in parallel on each BAM file
    coverage_results = RUN_COVERAGE(bam_files)

    // Merge all coverage JSON files into a single JSON output
    MERGE_COVERAGE(coverage_results.collect())
}