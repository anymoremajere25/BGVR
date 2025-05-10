nextflow.enable.dsl=2

params.sample_list = 'samples.txt'
params.region     = 'chr1:1-32'
params.mock = true  // Set to true to use mock commands for testing

workflow {
    // Define the channel for BAM files
    bamChannel = Channel.fromPath(params.sample_list)
                        .splitText()
                        .map { it.trim() }
    
    // Run the QC collection process
    qc_results = collectQC(bamChannel)
    
    // Merge QC results
    merged_qc = mergeQC(qc_results.collect())
    
    // Recalibrate each BAM with the merged QC data
    // Use combine to pair merged QC with each BAM file
    recalibrate(merged_qc, bamChannel)
}

process collectQC {
    tag "${bam_file}"
    
    input:
    val bam_file
    
    output:
    path "qc_${bam_file}.txt"
    
    script:
    if (params.mock) {
        """
        # Mock QC collection
        echo "Creating mock QC data for ${bam_file}"
        echo "bam_file: ${bam_file}" > qc_${bam_file}.txt
        echo "region: ${params.region}" >> qc_${bam_file}.txt
        echo "total_reads: 1000" >> qc_${bam_file}.txt
        echo "mapped_reads: 950" >> qc_${bam_file}.txt
        echo "coverage: 30.5" >> qc_${bam_file}.txt
        """
    } else {
        """
        /mnt/c/Users/trian/BGVR/chapter_06/experiment_66/target/debug/coverage_tool \
          --bam ${bam_file} \
          --region ${params.region} > qc_${bam_file}.txt
        """
    }
}

process mergeQC {
    input:
    path qc_files
    
    output:
    path "merged_qc.json"
    
    script:
    if (params.mock) {
        """
        # Mock QC merging
        echo "Merging QC files: ${qc_files}"
        echo "{" > merged_qc.json
        echo "  \"samples\": [" >> merged_qc.json
        
        # Process each QC file
        FIRST=true
        for file in ${qc_files}; do
            if [ "\$FIRST" = "true" ]; then
                FIRST=false
            else
                echo "," >> merged_qc.json
            fi
            
            # Extract sample name from filename
            SAMPLE=\$(basename \$file | sed 's/qc_\\(.*\\)\\.txt/\\1/')
            
            # Add sample data to JSON
            echo "    {" >> merged_qc.json
            echo "      \"sample\": \"\$SAMPLE\"," >> merged_qc.json
            echo "      \"region\": \"${params.region}\"," >> merged_qc.json
            echo "      \"total_reads\": 1000," >> merged_qc.json
            echo "      \"mapped_reads\": 950," >> merged_qc.json
            echo "      \"coverage\": 30.5" >> merged_qc.json
            echo "    }" >> merged_qc.json
        done
        
        echo "  ]" >> merged_qc.json
        echo "}" >> merged_qc.json
        """
    } else {
        """
        rust_merge_qc ${qc_files.join(" ")} merged_qc.json
        """
    }
}

process recalibrate {
    tag "${bam_file}"
    
    input:
    path merged_qc
    val bam_file
    
    output:
    path "recalibrated_${bam_file}"
    
    script:
    if (params.mock) {
        """
        # Mock recalibration
        echo "Recalibrating ${bam_file} using ${merged_qc}"
        echo "This is a mock recalibrated BAM file for ${bam_file}" > recalibrated_${bam_file}
        echo "Created using QC data from ${merged_qc}" >> recalibrated_${bam_file}
        echo "Original BAM: ${bam_file}" >> recalibrated_${bam_file}
        echo "Region: ${params.region}" >> recalibrated_${bam_file}
        """
    } else {
        """
        rust_recalibrator \
          --bam ${bam_file} \
          --qc ${merged_qc} \
          --out recalibrated_${bam_file}
        """
    }
}