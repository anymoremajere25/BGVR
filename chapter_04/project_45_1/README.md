## 4.6 TRANSCRIPTOMICS AND ALTERNATIVE SPLICING ALGORITHMS
### Explanation of the Rust Program for Splicing Graph Construction [project_45_1]

This Rust program constructs a splicing graph from genomic alignment data using parallel processing with Rayon.
### 1. Purpose of the Program

    Reads sequencing alignments (chromosomal positions of RNA reads).
    Builds a splicing graph, representing exon connectivity.
    Uses parallel processing (Rayon) to optimize graph construction.
    Saves the final splicing graph to a file.

### 2. Understanding the Code
(a) Data Structures
1️⃣ ExonSegment: Represents an exon

#[derive(Debug, Clone)]
struct ExonSegment {
    _start: usize,
    _end: usize,
}

    Represents an exon (a segment of a gene).
    _start and _end define its boundaries.
    The underscore _ avoids warnings if the field is unused.

2️⃣ Alignment: Represents a sequencing read

#[derive(Debug, Clone)]
struct Alignment {
    _chrom: String, // Chromosome name
    start: usize,   // Start position of alignment
    end: usize,     // End position of alignment
}

    Stores information about a read alignment.
    start and end determine the exon’s position.

3️⃣ SplicingGraph: Stores exon connections

#[derive(Debug, Default, Clone)]
struct SplicingGraph {
    adjacency: HashMap<(usize, usize), Vec<(ExonSegment, u64)>>,
}

    Represents a splicing graph as an adjacency list.
    adjacency is a HashMap where:
        Key: (start, end) exon position.
        Value: List of (ExonSegment, coverage), representing connected exons.

(b) Graph Operations
1️⃣ Adding Exon Junctions

fn add_junction(&mut self, exon_key: (usize, usize), target_exon: ExonSegment) {
    let coverage = 1;
    self.adjacency
        .entry(exon_key)
        .or_default()
        .push((target_exon, coverage));
}

    Adds a connection between exons.
    Uses entry().or_default() to ensure the key exists before adding.

2️⃣ Merging Graphs

fn merge(&mut self, other: SplicingGraph) {
    for (key, edges) in other.adjacency {
        self.adjacency.entry(key).or_default().extend(edges);
    }
}

    Merges another splicing graph into the current one.
    Uses entry().or_default() to combine adjacency lists.

(c) Processing Alignments
1️⃣ Processing a Batch of Alignments

fn process_alignment_chunk(batch: &[Alignment]) -> SplicingGraph {
    let mut local_graph = SplicingGraph::default();
    for align in batch {
        let exon_key = (align.start, align.end);
        let target_exon = ExonSegment {
            _start: align.start + 50,
            _end: align.end + 100,
        };
        local_graph.add_junction(exon_key, target_exon);
    }
    local_graph
}

    Processes a batch of alignments.
    Creates exon connections by shifting positions (+50, +100).

### 3. Parallel Graph Construction

let final_graph = chunks
    .into_par_iter()
    .map(|batch| process_alignment_chunk(&batch))
    .reduce(
        SplicingGraph::default,
        |mut acc, local_graph| {
            acc.merge(local_graph);
            acc
        },
    );

    Splits alignments into chunks.
    Processes chunks in parallel using Rayon.
    Merges all local graphs into a final splicing graph.

### 4. Writing the Graph to File

let out_file = File::create("partial_splicing_graph.bin")?;
let mut writer = BufWriter::new(out_file);
writer.write_all(b"Serialized splicing graph example\n")?;
writer.write_all(format!("{:#?}", final_graph).as_bytes())?;

    Writes the splicing graph to a binary file.

### 5. Example Output
(a) File Output (partial_splicing_graph.bin)

Serialized splicing graph example
SplicingGraph {
    adjacency: {
        (100, 200): [
            (ExonSegment { _start: 150, _end: 300 }, 1),
        ],
        (150, 300): [
            (ExonSegment { _start: 200, _end: 400 }, 1),
        ],
        (500, 700): [
            (ExonSegment { _start: 550, _end: 800 }, 1),
        ],
    },
}

    Shows the splicing graph structure.

(b) Console Output

Splicing graph has been successfully written.

    Confirms the file was saved.

### 6. Summary

✔ Builds a splicing graph from alignments.
✔ Uses parallel processing (Rayon).
✔ Efficiently merges graphs.
✔ Writes the graph to disk.
