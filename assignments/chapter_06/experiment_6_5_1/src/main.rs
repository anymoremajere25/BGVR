use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use env_logger;
use log::info;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::Arc;

#[derive(Clone, Debug)]
struct Interval {
    start: u64,
    end: u64,
    coverage: u64,
}

#[derive(Debug)]
struct IntervalTree {
    center: u64,
    intervals: Vec<Interval>,
    left: Option<Box<IntervalTree>>,
    right: Option<Box<IntervalTree>>,
}

#[derive(Parser, Debug)]
#[command(name = "genomic_interval_tree", about = "Tool for genomic coverage computation and queries")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Coverage {
        #[arg(long)]
        bam: String,
        #[arg(long)]
        out: String,
        #[arg(long, default_value_t = 50000)]
        chunk_size: u64,
    },
    Query {
        #[arg(long)]
        interval_file: String,
        #[arg(long)]
        query: String,
        #[arg(long)]
        output: String,
    },
}

impl IntervalTree {
    fn build(intervals: &[Interval]) -> Option<Box<IntervalTree>> {
        if intervals.is_empty() {
            return None;
        }
        let mut sorted = intervals.to_vec();
        sorted.sort_by_key(|iv| iv.start);
        let median_idx = sorted.len() / 2;
        let center = sorted[median_idx].start;
        let (left_intervals, mut right_candidates): (Vec<_>, Vec<_>) =
            sorted.into_iter().partition(|iv| iv.end < center);
        let center_intervals: Vec<Interval> = right_candidates
            .iter()
            .cloned()
            .filter(|iv| iv.start <= center && iv.end >= center)
            .collect();
        let right_intervals: Vec<Interval> = right_candidates
            .drain(..)
            .filter(|iv| iv.start > center)
            .collect();
        let left_tree = IntervalTree::build(&left_intervals);
        let right_tree = IntervalTree::build(&right_intervals);
        Some(Box::new(IntervalTree {
            center,
            intervals: center_intervals,
            left: left_tree,
            right: right_tree,
        }))
    }

    fn query(&self, qstart: u64, qend: u64) -> Vec<&Interval> {
        let mut results = Vec::new();
        for iv in &self.intervals {
            if iv.end >= qstart && iv.start <= qend {
                results.push(iv);
            }
        }
        if let Some(ref left_tree) = self.left {
            if qstart <= self.center {
                results.extend(left_tree.query(qstart, qend));
            }
        }
        if let Some(ref right_tree) = self.right {
            if qend >= self.center {
                results.extend(right_tree.query(qstart, qend));
            }
        }
        results
    }
}

fn parse_interval_tsv(line: &str) -> Result<Interval> {
    let parts: Vec<&str> = line.trim().split('\t').collect();
    if parts.len() != 3 {
        return Err(anyhow::anyhow!("Invalid TSV format, expected start\tend\tcoverage"));
    }
    let start: u64 = parts[0].parse().with_context(|| format!("Failed to parse start: {}", parts[0]))?;
    let end: u64 = parts[1].parse().with_context(|| format!("Failed to parse end: {}", parts[1]))?;
    let coverage: u64 = parts[2].parse().with_context(|| format!("Failed to parse coverage: {}", parts[2]))?;
    if start > end {
        return Err(anyhow::anyhow!("Start cannot exceed end"));
    }
    Ok(Interval { start, end, coverage })
}

fn parse_bam_line(line: &str) -> Result<(u64, u64)> {
    let parts: Vec<&str> = line.trim().split_whitespace().collect();
    if parts.len() != 2 {
        return Err(anyhow::anyhow!("Invalid BAM-like format, expected start end"));
    }
    let start: u64 = parts[0].parse().with_context(|| format!("Failed to parse start: {}", parts[0]))?;
    let end: u64 = parts[1].parse().with_context(|| format!("Failed to parse end: {}", parts[1]))?;
    if start > end {
        return Err(anyhow::anyhow!("Start cannot exceed end"));
    }
    Ok((start, end))
}

fn parse_query(s: &str) -> Result<(u64, u64)> {
    let range = s.split(':').last().ok_or_else(|| anyhow::anyhow!("Invalid query format"))?;
    let parts: Vec<&str> = range.split('-').collect();
    if parts.len() != 2 {
        return Err(anyhow::anyhow!("Invalid query format, expected start-end"));
    }
    let start: u64 = parts[0].parse().with_context(|| format!("Failed to parse start: {}", parts[0]))?;
    let end: u64 = parts[1].parse().with_context(|| format!("Failed to parse end: {}", parts[1]))?;
    if start > end {
        return Err(anyhow::anyhow!("Query start cannot exceed end"));
    }
    Ok((start, end))
}

fn compute_coverage(bam_path: &str, out_path: &str, chunk_size: u64) -> Result<()> {
    let file = File::open(bam_path).with_context(|| format!("Failed to open BAM-like file: {}", bam_path))?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    let mut current_start = None;
    let mut current_end = 0;
    let mut coverage = 0;

    for line in reader.lines() {
        let (start, end) = parse_bam_line(&line?)?;
        if current_start.is_none() {
            current_start = Some(start);
            current_end = end;
            coverage = 1;
        } else if start <= current_end + chunk_size && end > current_end {
            current_end = end;
            coverage += 1;
        } else {
            intervals.push(Interval {
                start: current_start.unwrap(),
                end: current_end,
                coverage,
            });
            current_start = Some(start);
            current_end = end;
            coverage = 1;
        }
    }
    if let Some(start) = current_start {
        intervals.push(Interval { start, end: current_end, coverage });
    }

    let mut output = File::create(out_path).with_context(|| format!("Failed to create output: {}", out_path))?;
    for iv in intervals {
        writeln!(output, "{}\t{}\t{}", iv.start, iv.end, iv.coverage)?;
    }
    Ok(())
}

fn query_intervals(interval_file: &str, query: &str, output: &str) -> Result<()> {
    let file = File::open(interval_file).with_context(|| format!("Failed to open interval file: {}", interval_file))?;
    let reader = BufReader::new(file);
    let intervals: Vec<Interval> = reader
        .lines()
        .map(|line| parse_interval_tsv(&line?))
        .collect::<Result<Vec<Interval>>>()?;

    let tree = IntervalTree::build(&intervals).ok_or_else(|| anyhow::anyhow!("No intervals provided"))?;
    info!("Interval tree built with {} intervals.", intervals.len());

    let (qstart, qend) = parse_query(query)?;
    let tree_arc = Arc::new(tree);
    let hits = tree_arc.query(qstart, qend);

    let mut output_file = File::create(output).with_context(|| format!("Failed to create output: {}", output))?;
    for iv in &hits {
        writeln!(output_file, "{}\t{}\t{}", iv.start, iv.end, iv.coverage)?;
    }
    info!("Query [{}, {}] => found {} intervals", qstart, qend, hits.len());
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Coverage { bam, out, chunk_size } => {
            compute_coverage(&bam, &out, chunk_size)?;
            info!("Computed coverage for {} to {}", bam, out);
        }
        Commands::Query { interval_file, query, output } => {
            query_intervals(&interval_file, &query, &output)?;
        }
    }
    Ok(())
}
