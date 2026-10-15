//! Catalogued operations of the `sequence` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("sequence");

pub(super) const ROWS: &[Capability] = &[
    D.row("align", "sequence, sequence", "Alignment", "cpu"),
    D.row("parse_fasta", "text", "list", "cpu"),
    D.row("write_fasta", "records", "text", "cpu"),
    D.row("kmer_counts", "sequence", "list", "cpu"),
    D.row(
        "msa",
        "sequences,scoring,refinement_passes,memory_limit",
        "list[str]",
        "cpu,memory",
    ),
    D.row("neighbor_joining", "labels,distances", "Tree", "cpu,memory"),
    D.row("upgma", "labels,distances", "Tree", "cpu,memory"),
];
