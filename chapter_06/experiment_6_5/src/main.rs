use anyhow::{Context, Result};
use clap::Parser;
use env_logger;
use log::{info};
use rayon::prelude::*;
use std::sync::{Arc};

/// An interval with start, end, and associated coverage.
#[derive(Clone, Debug)]
struct Interval {
    start: u64,
    end: u64,
    coverage: u64,
}

/// A node of the interval tree, containing a center point, intervals that overlap this center,
/// and optional left and right subtrees.
#[derive(Debug)]
struct IntervalTree {
    center: u64,
    intervals: Vec<Interval>,
    left: Option<Box<IntervalTree>>,
    right: Option<Box<IntervalTree>>,
}

/// Builds and queries an interval tree for coverage. Demonstrates parallel querying using rayon.
#[derive(Parser, Debug)]
#[command(name = "interval_query_tool", about = "A robust interval tree tool for genomic coverage queries")]
struct Cli {
    /// A list of intervals specified as "start-end:coverage" (for demonstration)
    #[arg(long, value_name = "INTERVALS", num_args = 1..)]
    intervals: Vec<String>,

    /// A list of query ranges specified as "start-end" (for demonstration)
    #[arg(long, value_name = "QUERIES", num_args = 1..)]
    queries: Vec<String>,
}

impl IntervalTree {
    /// Builds an interval tree from a slice of Interval structs. If no intervals exist, returns None.
    fn build(intervals: &[Interval]) -> Option<Box<IntervalTree>> {
        if intervals.is_empty() {
            return None;
        }

        // Sort intervals by their start coordinate
        let mut sorted = intervals.to_vec();
        sorted.sort_by_key(|iv| iv.start);

        // Choose a center around the median interval's start
        let median_idx = sorted.len() / 2;
        let center = sorted[median_idx].start;

        // Partition intervals into those strictly to the left vs. those that could be center or right
        let (left_intervals, mut right_candidates): (Vec<_>, Vec<_>) =
            sorted.into_iter().partition(|iv| iv.end < center);

        // Extract intervals that actually span the center
        let center_intervals: Vec<Interval> = right_candidates
            .iter()
            .cloned()
            .filter(|iv| iv.start <= center && iv.end >= center)
            .collect();

        // Remaining intervals go to the right subtree
        let right_intervals: Vec<Interval> = right_candidates
            .drain(..)
            .filter(|iv| iv.start > center)
            .collect();

        // Recursively build left and right subtrees
        let left_tree = IntervalTree::build(&left_intervals);
        let right_tree = IntervalTree::build(&right_intervals);

        Some(Box::new(IntervalTree {
            center,
            intervals: center_intervals,
            left: left_tree,
            right: right_tree,
        }))
    }

    /// Query the tree for all intervals overlapping the query range [qstart, qend].
    fn query(&self, qstart: u64, qend: u64) -> Vec<&Interval> {
        let mut results = Vec::new();

        // Collect intervals in this node that overlap [qstart, qend].
        for iv in &self.intervals {
            if iv.end >= qstart && iv.start <= qend {
                results.push(iv);
            }
        }

        // Traverse left subtree if the query range extends to or beyond this node’s center.
        if let Some(ref left_tree) = self.left {
            if qstart <= self.center {
                results.extend(left_tree.query(qstart, qend));
            }
        }

        // Traverse right subtree if the query range extends to or beyond this node’s center.
        if let Some(ref right_tree) = self.right {
            if qend >= self.center {
                results.extend(right_tree.query(qstart, qend));
            }
        }
        results
    }
}

/// Parses an interval string of the form "start-end:coverage" into an Interval struct.
/// Returns an error if parsing fails, making it easier to integrate with robust error handling.
fn parse_interval(s: &str) -> Result<Interval> {
    // Example: "10-20:30"
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return Err(anyhow::anyhow!("Invalid interval format, expected start-end:coverage"));
    }

    let coverage: u64 = parts[1]
        .parse()
        .with_context(|| format!("Failed to parse coverage in interval: {s}"))?;

    let range_parts: Vec<&str> = parts[0].split('-').collect();
    if range_parts.len() != 2 {
        return Err(anyhow::anyhow!("Invalid range format, expected start-end"));
    }

    let start: u64 = range_parts[0]
        .parse()
        .with_context(|| format!("Failed to parse start coordinate in interval: {s}"))?;
    let end: u64 = range_parts[1]
        .parse()
        .with_context(|| format!("Failed to parse end coordinate in interval: {s}"))?;

    if start > end {
        return Err(anyhow::anyhow!("Start coordinate cannot exceed end coordinate."));
    }

    Ok(Interval { start, end, coverage })
}

/// Parses a query string of the form "start-end" into (start, end).
/// Returns an error if parsing fails.
fn parse_query(s: &str) -> Result<(u64, u64)> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() != 2 {
        return Err(anyhow::anyhow!("Invalid query format, expected start-end"));
    }

    let start: u64 = parts[0]
        .parse()
        .with_context(|| format!("Failed to parse start coordinate in query: {s}"))?;
    let end: u64 = parts[1]
        .parse()
        .with_context(|| format!("Failed to parse end coordinate in query: {s}"))?;

    if start > end {
        return Err(anyhow::anyhow!("Query start cannot exceed query end."));
    }

    Ok((start, end))
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    // Parse input intervals from CLI arguments
    let intervals: Vec<Interval> = cli
        .intervals
        .iter()
        .map(|s| parse_interval(s))
        .collect::<Result<Vec<Interval>>>()?;

    // Build the interval tree. This tree is immutable after construction, making it safe to share.
    let tree = IntervalTree::build(&intervals).ok_or_else(|| anyhow::anyhow!("No intervals provided"))?;
    // Replaced info! with println!
    println!("Interval tree built with {} intervals.", intervals.len());

    // Parse query ranges from CLI arguments
    let queries: Vec<(u64, u64)> = cli
        .queries
        .iter()
        .map(|s| parse_query(s))
        .collect::<Result<Vec<(u64, u64)>>>()?;

    // Wrap the tree in an Arc to allow parallel access without data races
    let tree_arc = Arc::new(tree);

    // Perform queries in parallel using rayon
    queries.par_iter().for_each(|(qstart, qend)| {
        let hits = tree_arc.query(*qstart, *qend);
        // Replaced info! with println!
        println!("Query [{}, {}] => found {} intervals", qstart, qend, hits.len());
    });

    Ok(())
}