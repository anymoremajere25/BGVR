#!/usr/bin/env nextflow

/*
 * Epigenomic HMM Segmentation Pipeline
 * 
 * This pipeline generates synthetic epigenomic data and performs
 * HMM-based chromatin state segmentation using Rust
 */

nextflow.enable.dsl=2

// Default parameters
params.n_positions = 5000
params.n_states = 3
params.min_segment_length = 50
params.chromosome = "chr1"
params.start_position = 1000000
params.noise_level = 0.1
params.seed = 42
params.add_noise = false
params.visualize = true
params.max_iterations = 100
params.tolerance = 1e-6
params.parallel = false
params.outdir = "results"

// Input validation
if (params.n_states < 2) {
    error "Number of states must be at least 2"
}

if (params.n_positions < 100) {
    error "Number of positions must be at least 100"
}

log.info """\
    EPIGENOMIC HMM SEGMENTATION PIPELINE
    ===================================
    n_positions      : ${params.n_positions}
    n_states         : ${params.n_states}
    min_segment_length: ${params.min_segment_length}
    chromosome       : ${params.chromosome}
    start_position   : ${params.start_position}
    noise_level      : ${params.noise_level}
    seed             : ${params.seed}
    add_noise        : ${params.add_noise}
    visualize        : ${params.visualize}
    max_iterations   : ${params.max_iterations}
    tolerance        : ${params.tolerance}
    parallel         : ${params.parallel}
    outdir           : ${params.outdir}
    """

/*
 * Process 1: Generate synthetic epigenomic data
 */
process GENERATE_DATA {
    publishDir "${params.outdir}/data", mode: 'copy'
    
    input:
    val n_positions
    val n_states
    val min_segment_length
    val chromosome
    val start_position
    val noise_level
    val seed
    val add_noise
    val visualize
    
    output:
    path "epigenomic_data.csv", emit: data
    path "epigenomic_data_with_true_states.csv", emit: data_with_truth
    path "epigenomic_data_visualization.png", emit: plot, optional: true
    path "data_generation.log", emit: log
    
    script:
    def noise_flag = add_noise ? "--add-noise" : ""
    def viz_flag = visualize ? "--visualize" : ""
    
    """
    python3 generate_data.py \\
        --n-positions ${n_positions} \\
        --n-states ${n_states} \\
        --min-segment-length ${min_segment_length} \\
        --chromosome ${chromosome} \\
        --start-position ${start_position} \\
        --noise-level ${noise_level} \\
        --seed ${seed} \\
        --output epigenomic_data.csv \\
        ${noise_flag} \\
        ${viz_flag} \\
        2>&1 | tee data_generation.log
    """
}

/*
 * Process 2: Build Rust HMM segmentation tool
 */
process BUILD_RUST_TOOL {
    
    output:
    path "target/release/main", emit: executable
    path "build.log", emit: log
    
    script:
    """
    # Build the Rust project
    cargo build --release 2>&1 | tee build.log
    
    # Verify the executable was created
    if [ ! -f target/release/main ]; then
        echo "ERROR: Failed to build Rust executable" >&2
        exit 1
    fi
    
    echo "Successfully built Rust HMM segmentation tool"
    """
}

/*
 * Process 3: Run HMM segmentation
 */
process RUN_HMM_SEGMENTATION {
    publishDir "${params.outdir}/segmentation", mode: 'copy'
    
    input:
    path data_file
    path executable
    val n_states
    val max_iterations
    val tolerance
    val parallel
    
    output:
    path "segmentation_results.json", emit: results
    path "segmentation.log", emit: log
    
    script:
    def parallel_flag = parallel ? "--parallel" : ""
    
    """
    # Make executable runnable
    chmod +x ${executable}
    
    # Run HMM segmentation
    ./${executable} \\
        --input ${data_file} \\
        --output segmentation_results.json \\
        --states ${n_states} \\
        --max-iterations ${max_iterations} \\
        --tolerance ${tolerance} \\
        ${parallel_flag} \\
        2>&1 | tee segmentation.log
    """
}

/*
 * Process 4: Evaluate segmentation results
 */
process EVALUATE_RESULTS {
    publishDir "${params.outdir}/evaluation", mode: 'copy'
    
    input:
    path results_file
    path data_with_truth
    
    output:
    path "evaluation_report.html", emit: report
    path "evaluation_metrics.json", emit: metrics
    path "comparison_plots.png", emit: plots
    path "evaluation.log", emit: log
    
    script:
    """
    python3 evaluate_segmentation.py \\
        --results ${results_file} \\
        --truth ${data_with_truth} \\
        --output-report evaluation_report.html \\
        --output-metrics evaluation_metrics.json \\
        --output-plots comparison_plots.png \\
        2>&1 | tee evaluation.log
    """
}

/*
 * Process 5: Generate final report
 */
process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path data_log
    path build_log
    path segmentation_log
    path evaluation_log
    path evaluation_metrics
    path segmentation_results
    
    output:
    path "pipeline_report.html", emit: report
    path "pipeline_summary.json", emit: summary
    
    script:
    """
    python3 generate_pipeline_report.py \\
        --data-log ${data_log} \\
        --build-log ${build_log} \\
        --segmentation-log ${segmentation_log} \\
        --evaluation-log ${evaluation_log} \\
        --evaluation-metrics ${evaluation_metrics} \\
        --segmentation-results ${segmentation_results} \\
        --output-report pipeline_report.html \\
        --output-summary pipeline_summary.json \\
        --n-positions ${params.n_positions} \\
        --n-states ${params.n_states} \\
        --max-iterations ${params.max_iterations}
    """
}

/*
 * Main workflow
 */
workflow {
    // Generate synthetic data
    GENERATE_DATA(
        params.n_positions,
        params.n_states,
        params.min_segment_length,
        params.chromosome,
        params.start_position,
        params.noise_level,
        params.seed,
        params.add_noise,
        params.visualize
    )
    
    // Build Rust tool
    BUILD_RUST_TOOL()
    
    // Run HMM segmentation
    RUN_HMM_SEGMENTATION(
        GENERATE_DATA.out.data,
        BUILD_RUST_TOOL.out.executable,
        params.n_states,
        params.max_iterations,
        params.tolerance,
        params.parallel
    )
    
    // Evaluate results
    EVALUATE_RESULTS(
        RUN_HMM_SEGMENTATION.out.results,
        GENERATE_DATA.out.data_with_truth
    )
    
    // Generate final report
    GENERATE_REPORT(
        GENERATE_DATA.out.log,
        BUILD_RUST_TOOL.out.log,
        RUN_HMM_SEGMENTATION.out.log,
        EVALUATE_RESULTS.out.log,
        EVALUATE_RESULTS.out.metrics,
        RUN_HMM_SEGMENTATION.out.results
    )
}

/*
 * Workflow completion
 */
workflow.onComplete {
    log.info ( workflow.success ? """
        Pipeline completed successfully!
        Results are available in: ${params.outdir}
        
        Key outputs:
        - Data: ${params.outdir}/data/
        - Segmentation: ${params.outdir}/segmentation/
        - Evaluation: ${params.outdir}/evaluation/
        - Final Report: ${params.outdir}/pipeline_report.html
        """ : """
        Pipeline failed!
        Check the logs for error details.
        """ )
}
