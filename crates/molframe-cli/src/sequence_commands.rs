//! Standalone sequence file workflows backed by `molframe-seq`.

use crate::exit::Exit;
use crate::report::Context;
use crate::{
    KmerOperation, MatrixChoice, PairwiseMode, SequenceCommand, SequenceFormat, TreeMethod,
};
use molframe::seq::{Alignment, FastaRecord};
use std::io::Read as _;

pub(crate) fn execute(command: SequenceCommand, context: Context) -> Exit {
    let result = match command {
        SequenceCommand::Extract { input, chemistry } => {
            return crate::chemistry::with_ccd(
                chemistry.ccd.as_deref(),
                chemistry.ccd_version.as_deref(),
                context,
                |ccd, version, context| crate::inspect::sequence(&input, ccd, version, context),
            );
        }
        SequenceCommand::Convert { input, from, to } => convert(&input, from, to),
        SequenceCommand::Align {
            input,
            mode,
            matrix,
            gap_open,
            gap_extend,
        } => align(&input, mode, matrix, gap_open, gap_extend),
        SequenceCommand::Msa {
            input,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            refinement_passes,
        } => msa(
            &input,
            molframe::seq::Scoring {
                match_score,
                mismatch_score,
                gap_open,
                gap_extend,
            },
            refinement_passes,
        ),
        SequenceCommand::Kmer {
            input,
            k,
            operation,
            window,
        } => kmer(&input, k, operation, window),
        SequenceCommand::Tree { input, method } => tree(&input, method),
    };
    match result {
        Ok(text) => {
            context.result(&text);
            Exit::Success
        }
        Err(error) => {
            eprintln!("{error}");
            Exit::Failure
        }
    }
}

fn convert(input: &str, from: SequenceFormat, to: SequenceFormat) -> Result<String, String> {
    let text = read_text(input)?;
    let document = molframe::seq::read_sequence(&text, sequence_format(from)).map_err(display)?;
    molframe::seq::write_sequence(&document, sequence_format(to)).map_err(display)
}

fn align(
    input: &str,
    mode: PairwiseMode,
    matrix: MatrixChoice,
    gap_open: i32,
    gap_extend: i32,
) -> Result<String, String> {
    let text = read_text(input)?;
    let records = fasta_records(&text)?;
    if records.len() != 2 {
        return Err("pairwise alignment requires exactly two FASTA records".to_owned());
    }
    let matrix = substitution_matrix(matrix)?;
    let left = &records[0].sequence;
    let right = &records[1].sequence;
    let alignment = match mode {
        PairwiseMode::Global => {
            molframe::seq::global_matrix(left, right, &matrix, gap_open, gap_extend)
        }
        PairwiseMode::Local => {
            molframe::seq::local_matrix(left, right, &matrix, gap_open, gap_extend)
        }
        PairwiseMode::SemiGlobal => {
            molframe::seq::semi_global_matrix(left, right, &matrix, gap_open, gap_extend)
        }
    }
    .map_err(display)?;
    let (left_aligned, right_aligned) = aligned_sequences(left, right, &alignment)?;
    let output = vec![
        FastaRecord {
            id: records[0].id.clone(),
            description: records[0].description.clone(),
            sequence: left_aligned,
        },
        FastaRecord {
            id: records[1].id.clone(),
            description: records[1].description.clone(),
            sequence: right_aligned,
        },
    ];
    Ok(molframe::seq::write_fasta(&output))
}

fn msa(input: &str, scoring: molframe::seq::Scoring, passes: usize) -> Result<String, String> {
    let text = read_text(input)?;
    let mut records = fasta_records(&text)?;
    if records.is_empty() {
        return Err("multiple-sequence alignment requires FASTA records".to_owned());
    }
    let sequences: Vec<&[u8]> = records
        .iter()
        .map(|record| record.sequence.as_slice())
        .collect();
    let options = molframe::seq::MsaOptions::progressive(scoring).with_refinement_passes(passes);
    let aligned = molframe::seq::progressive_msa(&sequences, options).map_err(display)?;
    for (record, sequence) in records.iter_mut().zip(aligned) {
        record.sequence = sequence;
    }
    Ok(molframe::seq::write_fasta(&records))
}

fn kmer(
    input: &str,
    k: usize,
    operation: KmerOperation,
    window: Option<usize>,
) -> Result<String, String> {
    let text = read_text(input)?;
    let records = fasta_records(&text)?;
    let Some(record) = records.first() else {
        return Err("k-mer analysis requires at least one FASTA record".to_owned());
    };
    match operation {
        KmerOperation::Counts => Ok(molframe::seq::kmer_counts(&record.sequence, k)
            .into_iter()
            .map(|(word, count)| format!("{}\t{count}", String::from_utf8_lossy(&word)))
            .collect::<Vec<_>>()
            .join("\n")),
        KmerOperation::Minimizers => {
            let Some(window) = window else {
                return Err("minimizers require --window".to_owned());
            };
            Ok(molframe::seq::minimizers(&record.sequence, k, window)
                .into_iter()
                .map(|(position, word)| format!("{position}\t{}", String::from_utf8_lossy(&word)))
                .collect::<Vec<_>>()
                .join("\n"))
        }
    }
}

fn tree(input: &str, method: TreeMethod) -> Result<String, String> {
    let text = read_text(input)?;
    let (labels, distances) = parse_distance_matrix(&text)?;
    let label_refs: Vec<&str> = labels.iter().map(String::as_str).collect();
    let tree = match method {
        TreeMethod::NeighborJoining => molframe::seq::neighbor_joining(&label_refs, &distances),
        TreeMethod::Upgma => molframe::seq::upgma(&label_refs, &distances),
    };
    tree.map(|value| value.to_newick())
        .ok_or_else(|| "distance matrix cannot produce a tree".to_owned())
}

const fn sequence_format(format: SequenceFormat) -> molframe::seq::SequenceFormat {
    match format {
        SequenceFormat::Fasta => molframe::seq::SequenceFormat::Fasta,
        SequenceFormat::Fastq => molframe::seq::SequenceFormat::Fastq,
        SequenceFormat::A2m => molframe::seq::SequenceFormat::A2m,
        SequenceFormat::A3m => molframe::seq::SequenceFormat::A3m,
        SequenceFormat::Clustal => molframe::seq::SequenceFormat::Clustal,
        SequenceFormat::Stockholm => molframe::seq::SequenceFormat::Stockholm,
        SequenceFormat::Phylip => molframe::seq::SequenceFormat::Phylip,
    }
}

fn fasta_records(text: &str) -> Result<Vec<FastaRecord>, String> {
    match molframe::seq::read_sequence(text, molframe::seq::SequenceFormat::Fasta)
        .map_err(display)?
    {
        molframe::seq::SequenceDocument::Records(records) => Ok(records),
        _ => Err("FASTA dispatcher returned an incompatible document".to_owned()),
    }
}

fn aligned_sequences(
    left: &[u8],
    right: &[u8],
    alignment: &Alignment,
) -> Result<(Vec<u8>, Vec<u8>), String> {
    let mut left_aligned = Vec::with_capacity(alignment.columns.len());
    let mut right_aligned = Vec::with_capacity(alignment.columns.len());
    for column in &alignment.columns {
        left_aligned.push(project_column(left, column.left)?);
        right_aligned.push(project_column(right, column.right)?);
    }
    Ok((left_aligned, right_aligned))
}

fn project_column(sequence: &[u8], index: Option<usize>) -> Result<u8, String> {
    match index {
        Some(index) => sequence
            .get(index)
            .copied()
            .ok_or_else(|| "native alignment returned an invalid sequence position".to_owned()),
        None => Ok(b'-'),
    }
}

fn parse_distance_matrix(text: &str) -> Result<(Vec<String>, Vec<Vec<f64>>), String> {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let Some(header) = lines.next() else {
        return Err("distance matrix is empty".to_owned());
    };
    let labels: Vec<String> = header.split('\t').map(str::to_owned).collect();
    if labels.len() < 2 || labels.iter().any(String::is_empty) {
        return Err("distance matrix header requires at least two tab-separated labels".to_owned());
    }
    let distances = lines
        .map(|line| {
            line.split('\t')
                .map(|value| {
                    value
                        .parse::<f64>()
                        .map_err(|_| format!("invalid distance value: {value}"))
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    if distances.len() != labels.len() || distances.iter().any(|row| row.len() != labels.len()) {
        return Err("distance matrix must be square and match its header".to_owned());
    }
    Ok((labels, distances))
}

fn read_text(input: &str) -> Result<String, String> {
    if input == "-" {
        let mut text = String::new();
        std::io::stdin()
            .read_to_string(&mut text)
            .map_err(|error| format!("could not read standard input: {error}"))?;
        Ok(text)
    } else {
        std::fs::read_to_string(input).map_err(|error| format!("could not read {input}: {error}"))
    }
}

fn display(error: impl std::fmt::Display) -> String {
    error.to_string()
}

fn substitution_matrix(choice: MatrixChoice) -> Result<molframe::seq::SubstitutionMatrix, String> {
    match choice {
        MatrixChoice::Blosum62 => Ok(molframe::seq::blosum62()),
        MatrixChoice::Identity => {
            molframe::seq::load_matrix(molframe::seq::MatrixProfile::Identity).map_err(display)
        }
        MatrixChoice::Nuc44 => {
            molframe::seq::load_matrix(molframe::seq::MatrixProfile::Nuc44).map_err(display)
        }
    }
}

#[cfg(test)]
#[path = "sequence_commands_tests.rs"]
mod tests;
