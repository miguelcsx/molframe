//! Typed views over the free-text `REMARK` records that carry data.
//!
//! `REMARK 2` holds the resolution, `REMARK 350` the biological assemblies and
//! `REMARK 465` the residues missing from the model.

use super::PdbHeaders;
use super::typed::{ResidueId, residue_id};
use crate::fixed;

/// A residue the model leaves out, from `REMARK 465`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MissingResidue {
    /// Model number for a multi-model entry.
    pub model: Option<u32>,
    /// The missing residue.
    pub residue: ResidueId,
}

/// One rigid operation of a biological assembly, from `BIOMT1..3`.
#[derive(Clone, Debug, PartialEq)]
pub struct Biomt {
    /// Operator serial number.
    pub serial: u32,
    /// Row-major rotation.
    pub rotation: [[f64; 3]; 3],
    /// Translation in ångström.
    pub translation: [f64; 3],
}

/// Operations applied to a set of chains.
#[derive(Clone, Debug, PartialEq)]
pub struct BiomtGroup {
    /// Chains the operations apply to.
    pub chains: Vec<Box<str>>,
    /// Operations in serial order.
    pub operations: Vec<Biomt>,
}

/// One biological assembly, from `REMARK 350`.
#[derive(Clone, Debug, PartialEq)]
pub struct Biomolecule {
    /// Assembly number.
    pub id: u32,
    /// One group per `APPLY THE FOLLOWING TO CHAINS` block.
    pub groups: Vec<BiomtGroup>,
}

impl PdbHeaders {
    /// The text of every `REMARK <number>` line, after the number.
    fn remark_text<'a>(&'a self, number: &'a str) -> impl Iterator<Item = &'a str> {
        self.named("REMARK").filter_map(move |record| {
            let line = record.line();
            (fixed::text(line, 8, 10) == number).then(|| fixed::raw(line, 11, line.len()))
        })
    }

    /// The resolution in ångström that `REMARK 2` deposits, when it gives one.
    #[must_use]
    pub fn resolution(&self) -> Option<f64> {
        self.remark_text("2").find_map(|text| {
            let rest = text.trim().strip_prefix("RESOLUTION.")?;
            rest.split_whitespace().next()?.parse().ok()
        })
    }

    /// The residues `REMARK 465` lists as missing from the model.
    #[must_use]
    pub fn missing_residues(&self) -> Vec<MissingResidue> {
        let mut in_table = false;
        let mut missing = Vec::new();
        for record in self.named("REMARK") {
            let line = record.line();
            if fixed::text(line, 8, 10) != "465" {
                continue;
            }
            if line.contains("RES C SSSEQI") {
                in_table = true;
                continue;
            }
            if !in_table {
                continue;
            }
            // Single-model entries leave the model column empty.
            let model = fixed::integer(line, 12, 14).and_then(|n| u32::try_from(n).ok());
            if let Some(residue) = residue_id(line, (16, 18), 20, (22, 26), 27) {
                missing.push(MissingResidue { model, residue });
            }
        }
        missing
    }

    /// The biological assemblies `REMARK 350` defines.
    #[must_use]
    pub fn biomolecules(&self) -> Vec<Biomolecule> {
        let mut assemblies: Vec<Biomolecule> = Vec::new();
        let mut continuing = false;
        for text in self.remark_text("350") {
            let text = text.trim();
            if let Some(id) = text.strip_prefix("BIOMOLECULE:") {
                continuing = false;
                if let Ok(id) = id.trim().parse() {
                    assemblies.push(Biomolecule {
                        id,
                        groups: Vec::new(),
                    });
                }
            } else if let Some(chains) = chain_list(text) {
                let Some(assembly) = assemblies.last_mut() else {
                    continue;
                };
                // `AND CHAINS:` extends the group a long list began.
                let append = continuing && text.starts_with("AND");
                match assembly.groups.last_mut() {
                    Some(group) if append => group.chains.extend(chains),
                    _ => assembly.groups.push(BiomtGroup {
                        chains,
                        operations: Vec::new(),
                    }),
                }
                continuing = true;
            } else if let Some(row) = biomt_row(text) {
                continuing = false;
                if let Some(group) = assemblies.last_mut().and_then(|a| a.groups.last_mut()) {
                    push_row(group, row);
                }
            }
        }
        assemblies
    }
}

/// The chain names a `APPLY THE FOLLOWING TO CHAINS:` or `AND CHAINS:` line holds.
fn chain_list(text: &str) -> Option<Vec<Box<str>>> {
    let (_, list) = text.split_once("CHAINS:")?;
    if !(text.starts_with("APPLY THE FOLLOWING TO") || text.starts_with("AND")) {
        return None;
    }
    Some(
        list.split(',')
            .map(str::trim)
            .filter(|chain| !chain.is_empty())
            .map(Box::from)
            .collect(),
    )
}

/// One `BIOMTk` row: its index within the operator, serial, and four numbers.
fn biomt_row(text: &str) -> Option<(usize, u32, [f64; 4])> {
    let mut fields = text.split_whitespace();
    let tag = fields.next()?.strip_prefix("BIOMT")?;
    let row: usize = tag.parse().ok()?;
    let serial = fields.next()?.parse().ok()?;
    let mut numbers = [0.0; 4];
    for slot in &mut numbers {
        *slot = fields.next()?.parse().ok()?;
    }
    (1..=3).contains(&row).then_some((row - 1, serial, numbers))
}

fn push_row(group: &mut BiomtGroup, (row, serial, numbers): (usize, u32, [f64; 4])) {
    if group
        .operations
        .last()
        .is_none_or(|last| last.serial != serial)
    {
        group.operations.push(Biomt {
            serial,
            rotation: [[0.0; 3]; 3],
            translation: [0.0; 3],
        });
    }
    if let Some(operation) = group.operations.last_mut() {
        operation.rotation[row] = [numbers[0], numbers[1], numbers[2]];
        operation.translation[row] = numbers[3];
    }
}

#[cfg(test)]
#[path = "remark_tests.rs"]
mod tests;
