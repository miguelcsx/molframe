//! Declarative standalone sequence workflows.

use super::workflows::CcdArguments;
use clap::{Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Subcommand, Debug)]
pub(crate) enum SequenceCommand {
    /// Extract one-letter polymer sequences from a structure.
    Extract {
        input: PathBuf,
        #[command(flatten)]
        chemistry: CcdArguments,
    },
    /// Convert between sequence and alignment text formats.
    Convert {
        /// Input path, or `-` for standard input.
        input: String,
        #[arg(long, value_enum)]
        from: SequenceFormat,
        #[arg(long, value_enum)]
        to: SequenceFormat,
    },
    /// Align the first two records in a sequence file.
    Align {
        /// FASTA input path, or `-` for standard input.
        input: String,
        #[arg(long, value_enum, default_value_t = PairwiseMode::Global)]
        mode: PairwiseMode,
        #[arg(long, value_enum, default_value_t = MatrixChoice::Blosum62)]
        matrix: MatrixChoice,
        #[arg(long, default_value_t = -11)]
        gap_open: i32,
        #[arg(long, default_value_t = -1)]
        gap_extend: i32,
    },
    /// Build a deterministic progressive multiple-sequence alignment.
    Msa {
        /// FASTA input path, or `-` for standard input.
        input: String,
        #[arg(long, default_value_t = 1)]
        match_score: i32,
        #[arg(long, default_value_t = -1)]
        mismatch_score: i32,
        #[arg(long, default_value_t = -2)]
        gap_open: i32,
        #[arg(long, default_value_t = -1)]
        gap_extend: i32,
        #[arg(long, default_value_t = 0)]
        refinement_passes: usize,
    },
    /// Count exact k-mers or emit minimizers for the first FASTA record.
    Kmer {
        /// FASTA input path, or `-` for standard input.
        input: String,
        #[arg(long)]
        k: usize,
        #[arg(long, value_enum, default_value_t = KmerOperation::Counts)]
        operation: KmerOperation,
        /// Window length required for minimizers.
        #[arg(long, required_if_eq("operation", "minimizers"))]
        window: Option<usize>,
    },
    /// Build a tree from a labelled tab-separated distance matrix.
    Tree {
        /// Matrix path, or `-` for standard input; first row contains labels.
        input: String,
        #[arg(long, value_enum, default_value_t = TreeMethod::NeighborJoining)]
        method: TreeMethod,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(crate) enum SequenceFormat {
    Fasta,
    Fastq,
    A2m,
    A3m,
    Clustal,
    Stockholm,
    Phylip,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub(crate) enum PairwiseMode {
    #[default]
    Global,
    Local,
    SemiGlobal,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub(crate) enum MatrixChoice {
    #[default]
    Blosum62,
    Identity,
    Nuc44,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub(crate) enum KmerOperation {
    #[default]
    Counts,
    Minimizers,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub(crate) enum TreeMethod {
    #[default]
    NeighborJoining,
    Upgma,
}
