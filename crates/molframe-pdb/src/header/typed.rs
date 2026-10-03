//! Typed views over the preserved `SEQRES`, `SSBOND` and `LINK` records.
//!
//! The deposited lines stay the source of truth, so these views are computed on
//! demand and writing a structure back out reproduces the header unchanged.

use super::PdbHeaders;
use crate::fixed;

/// A residue named the way connectivity records name it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResidueId {
    /// Component name, trimmed.
    pub name: Box<str>,
    /// Chain identifier; empty for a file that leaves the column blank.
    pub chain: Box<str>,
    /// Residue sequence number.
    pub sequence: i32,
    /// Insertion code, when one is given.
    pub insertion: Option<char>,
}

/// The residue sequence one chain declares in its `SEQRES` records.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqRes {
    /// Chain identifier.
    pub chain: Box<str>,
    /// Residue count the first record declares.
    pub declared: Option<u32>,
    /// Component names in sequence order.
    pub residues: Vec<Box<str>>,
}

/// A disulfide bond declared by `SSBOND`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SsBond {
    /// First bonded residue.
    pub first: ResidueId,
    /// Second bonded residue.
    pub second: ResidueId,
}

/// A covalent link declared by `LINK`.
#[derive(Clone, Debug, PartialEq)]
pub struct Link {
    /// First atom name.
    pub first_atom: Box<str>,
    /// Residue of the first atom.
    pub first: ResidueId,
    /// Second atom name.
    pub second_atom: Box<str>,
    /// Residue of the second atom.
    pub second: ResidueId,
    /// Deposited link length in ångström, when given.
    pub distance: Option<f64>,
}

/// Reads a residue from fixed columns: name, chain, sequence, insertion code.
pub(super) fn residue_id(
    line: &str,
    name: (usize, usize),
    chain: usize,
    sequence: (usize, usize),
    insertion: usize,
) -> Option<ResidueId> {
    let sequence = i32::try_from(fixed::integer(line, sequence.0, sequence.1)?).ok()?;
    Some(ResidueId {
        name: fixed::text(line, name.0, name.1).into(),
        chain: fixed::text(line, chain, chain).into(),
        sequence,
        insertion: fixed::text(line, insertion, insertion).chars().next(),
    })
}

impl PdbHeaders {
    /// The `SEQRES` residue sequences, one per chain in first-appearance order.
    #[must_use]
    pub fn seqres(&self) -> Vec<SeqRes> {
        let mut chains: Vec<SeqRes> = Vec::new();
        for record in self.named("SEQRES") {
            let line = record.line();
            let chain = fixed::text(line, 12, 12);
            if !chains.iter().any(|seen| &*seen.chain == chain) {
                chains.push(SeqRes {
                    chain: chain.into(),
                    declared: fixed::integer(line, 14, 17).and_then(|n| u32::try_from(n).ok()),
                    residues: Vec::new(),
                });
            }
            let Some(entry) = chains.iter_mut().find(|seen| &*seen.chain == chain) else {
                continue;
            };
            // Up to thirteen three-column names, four columns apart.
            for slot in 0..13 {
                let from = 20 + 4 * slot;
                let name = fixed::text(line, from, from + 2);
                if !name.is_empty() {
                    entry.residues.push(name.into());
                }
            }
        }
        chains
    }

    /// The `SSBOND` disulfide bonds, in deposited order.
    #[must_use]
    pub fn ssbonds(&self) -> Vec<SsBond> {
        self.named("SSBOND")
            .filter_map(|record| {
                let line = record.line();
                Some(SsBond {
                    first: residue_id(line, (12, 14), 16, (18, 21), 22)?,
                    second: residue_id(line, (26, 28), 30, (32, 35), 36)?,
                })
            })
            .collect()
    }

    /// The `LINK` records, in deposited order.
    #[must_use]
    pub fn links(&self) -> Vec<Link> {
        self.named("LINK")
            .filter_map(|record| {
                let line = record.line();
                Some(Link {
                    first_atom: fixed::text(line, 13, 16).into(),
                    first: residue_id(line, (18, 20), 22, (23, 26), 27)?,
                    second_atom: fixed::text(line, 43, 46).into(),
                    second: residue_id(line, (48, 50), 52, (53, 56), 57)?,
                    distance: fixed::real(line, 74, 78),
                })
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "typed_tests.rs"]
mod tests;
