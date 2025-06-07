use clap::{Arg, Command};
use chipseq_peak_caller::{PeakCaller, PeakCallingConfig};
use std::process;

fn main() {
    let matches = Command::new("ChIP-Seq Peak Caller")
        .version("0.1.0")
        .author("Your Name <your.email@example.com>")
        .about("A simple peak-calling tool for ChIP-Seq data")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("Input coverage file in BedGraph format")
                .required(true)
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output peaks file in BED format")
                .default_value("peaks.bed")
        )
        .arg(
            Arg::new("fold-enrichment")
                .short('f')
                .long("fold-enrichment")
                .value_name("FLOAT")
                .help("Minimum fold enrichment threshold")
                .default_value("2.0")
        )
        .arg(
            Arg::new("background-window")
                .short('w')
                .long("background-window")
                .value_name("INT")
                .help("Background window size in base pairs")
                .default_value("10000")
        )
        .arg(
            Arg::new("min-peak-length")
                .short('l')
                .long("min-peak-length")
                .value_name("INT")
                .help("Minimum peak length in base pairs")
                .default_value("150")
        )
        .arg(
            Arg::new("merge-distance")
                .short('m')
                .long("merge-distance")
                .value_name("INT")
                .help("Merge peaks within this distance")
                .default_value("200")
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .help("Enable verbose output")
                .action(clap::ArgAction::SetTrue)
        )
        .get_matches();

    // Parse command line arguments
    let input_file = matches.get_one::<String>("input").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let verbose = matches.get_flag("verbose");
    
    let fold_enrichment: f64 = matches.get_one::<String>("fold-enrichment")
        .unwrap()
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Error: Invalid fold enrichment value");
            process::exit(1);
        });
    
    let background_window: u32 = matches.get_one::<String>("background-window")
        .unwrap()
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Error: Invalid background window size");
            process::exit(1);
        });
    
    let min_peak_length: u32 = matches.get_one::<String>("min-peak-length")
        .unwrap()
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Error: Invalid minimum peak length");
            process::exit(1);
        });
    
    let merge_distance: u32 = matches.get_one::<String>("merge-distance")
        .unwrap()
        .parse()
        .unwrap_or_else(|_| {
            eprintln!("Error: Invalid merge distance");
            process::exit(1);
        });

    // Create configuration
    let config = PeakCallingConfig {
        min_fold_enrichment: fold_enrichment,
        background_window,
        min_peak_length,
        merge_distance,
    };

    if verbose {
        println!("Configuration:");
        println!("  Input file: {}", input_file);
        println!("  Output file: {}", output_file);
        println!("  Fold enrichment threshold: {}", config.min_fold_enrichment);
        println!("  Background window: {} bp", config.background_window);
        println!("  Minimum peak length: {} bp", config.min_peak_length);
        println!("  Merge distance: {} bp", config.merge_distance);
        println!();
    }

    // Create peak caller and load data
    let mut caller = PeakCaller::with_config(config);
    
    if verbose {
        println!("Loading coverage data from {}...", input_file);
    }
    
    if let Err(e) = caller.load_coverage_from_bedgraph(input_file) {
        eprintln!("Error loading coverage data: {}", e);
        process::exit(1);
    }

    if verbose {
        println!("Calling peaks...");
    }
    
    // Call peaks
    let peaks = caller.call_peaks();
    
    if verbose {
        println!("Found {} peaks", peaks.len());
    }

    // Write output
    if let Err(e) = caller.write_peaks_to_bed(&peaks, output_file) {
        eprintln!("Error writing output: {}", e);
        process::exit(1);
    }

    // Print summary
    let summary = caller.get_summary(&peaks);
    if verbose {
        println!("\n{}", summary);
    } else {
        println!("Peak calling completed. Found {} peaks in {} bp total coverage.", 
                summary.total_peaks, summary.total_coverage);
    }
    
    println!("Results written to: {}", output_file);
}