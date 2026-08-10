//! The invariants a structure must satisfy.
//!
//! These are checked when a structure is built and whenever an edit is
//! committed, so a violation is caught at the one place that could have caused
//! it rather than surfacing later as an implausible number. Every violation is
//! reported rather than only the first, because a structure with a broken
//! hierarchy usually has more than one thing wrong with it.
//!
//! A violation here is a defect in whatever built the structure. No input file,
//! however malformed, should be able to produce one — a reader that cannot make
//! sense of a file emits a diagnostic and declines, rather than building
//! something inconsistent.

use super::data::{CoordinateStore, StructureData};
use crate::diagnostic::{Code, Diagnostic, Diagnostics};
use crate::index::ResidueIndex;

/// Checks every invariant, returning what was violated.
///
/// An empty result means the structure is sound.
#[must_use]
pub fn validate(data: &StructureData) -> Vec<Diagnostic> {
    let mut findings = Diagnostics::new();
    check_child_ranges(data, &mut findings);
    check_atom_coverage(data, &mut findings);
    check_entity_references(data, &mut findings);
    check_bonds(data, &mut findings);
    check_annotations(data, &mut findings);
    check_coordinate_counts(data, &mut findings);
    check_coordinates(data, &mut findings);
    check_occupancies(data, &mut findings);
    findings.finish()
}

fn check_annotations(data: &StructureData, findings: &mut Diagnostics) {
    let atoms = data.atom_count();
    for (name, column) in data.annotations.iter() {
        if column.len() != atoms {
            findings.push(
                Diagnostic::new(Code::E3011)
                    .with_context("annotation", name)
                    .with_context("column rows", column.len().to_string())
                    .with_context("atoms", atoms.to_string()),
            );
        }
    }
}

fn check_bonds(data: &StructureData, findings: &mut Diagnostics) {
    let atom_count = data.chunks.last().map_or(0, |chunk| chunk.atoms().end);
    for (position, bond) in data.bonds.iter().enumerate() {
        if bond.atom_a == bond.atom_b {
            findings.push(Diagnostic::new(Code::E3007).with_context("bond", position.to_string()));
        }
        for endpoint in [bond.atom_a, bond.atom_b] {
            if endpoint.get() >= atom_count {
                findings.push(
                    Diagnostic::new(Code::E3006)
                        .with_context("bond", position.to_string())
                        .with_context("atom", endpoint.to_string()),
                );
            }
        }
    }
}

fn check_coordinates(data: &StructureData, findings: &mut Diagnostics) {
    if let CoordinateStore::Ragged { models } = &data.coords {
        for (model, structure) in models.iter().enumerate() {
            for finding in validate(structure.data()) {
                findings.push(finding.with_context("ragged model", model.to_string()));
            }
        }
        return;
    }
    for model in 0..data.coords.model_count() {
        let Ok(model) = u32::try_from(model) else {
            findings.push(Diagnostic::new(Code::E3001).with_context("model", "exceeds u32"));
            return;
        };
        let Some(block) = data.coords.block(crate::index::ModelIndex::new(model)) else {
            continue;
        };
        for chunk in data.chunks.iter() {
            for local in 0..chunk.len() {
                if !chunk.has_position(local) {
                    continue;
                }
                let atom = chunk.atoms().start + local;
                if block
                    .as_slice()
                    .get(atom as usize)
                    .is_some_and(|position| position.iter().any(|value| !value.is_finite()))
                {
                    findings.push(
                        Diagnostic::new(Code::E3008)
                            .with_context("model", model.to_string())
                            .with_context("atom", atom.to_string()),
                    );
                }
            }
        }
    }
}

/// Every parent's children must lie inside the child table, and siblings must
/// not overlap or run backwards.
fn check_child_ranges(data: &StructureData, findings: &mut Diagnostics) {
    let topology = &data.topology;

    let mut expected_chain = 0u32;
    let shared_topology = matches!(&data.coords, CoordinateStore::Dense { .. });
    let mut shared_range = None;
    for model in topology.models.iter() {
        let Some(range) = topology.models.chains(model) else {
            continue;
        };
        if range.end as usize > topology.chains.len() {
            findings.push(
                Diagnostic::new(Code::E3002)
                    .with_message("a model's chains run past the chain table")
                    .with_context("model", model.to_string()),
            );
        }
        if shared_topology {
            if let Some(expected) = &shared_range {
                if expected != &range {
                    findings.push(
                        Diagnostic::new(Code::E3003)
                            .with_message("dense models do not share one chain range")
                            .with_context("model", model.to_string()),
                    );
                }
            } else {
                shared_range = Some(range.clone());
            }
        } else if range.start < expected_chain {
            findings.push(
                Diagnostic::new(Code::E3003)
                    .with_message("models share chains")
                    .with_context("model", model.to_string()),
            );
        }
        if !shared_topology {
            expected_chain = range.end;
        }
    }

    let mut expected_residue = 0u32;
    for chain in topology.chains.iter() {
        let Some(range) = topology.chains.residues(chain) else {
            continue;
        };
        if range.end as usize > topology.residues.len() {
            findings.push(
                Diagnostic::new(Code::E3002)
                    .with_message("a chain's residues run past the residue table")
                    .with_context("chain", chain.to_string()),
            );
        }
        if range.start < expected_residue {
            findings.push(
                Diagnostic::new(Code::E3003)
                    .with_message("chains share residues")
                    .with_context("chain", chain.to_string()),
            );
        }
        expected_residue = range.end;
    }
}

/// Every atom must belong to exactly one residue, which means the residues'
/// atom ranges must tile the atom order without gaps or overlaps.
fn check_atom_coverage(data: &StructureData, findings: &mut Diagnostics) {
    let residues = &data.topology.residues;
    let mut expected = 0u32;
    let Ok(residue_count) = u32::try_from(residues.len()) else {
        findings.push(Diagnostic::new(Code::E3001).with_context("residues", "exceeds u32"));
        return;
    };
    for position in 0..residue_count {
        let Some(range) = residues.atoms(ResidueIndex::new(position)) else {
            continue;
        };
        if range.start != expected {
            findings.push(
                Diagnostic::new(Code::E3004)
                    .with_message("residues do not cover the atoms contiguously")
                    .with_context("residue", position.to_string())
                    .with_context("expected first atom", expected.to_string()),
            );
        }
        expected = range.end.max(expected);
    }

    let atoms_in_chunks = data.chunks.last().map_or(0, |chunk| chunk.atoms().end);
    if !residues.is_empty() && expected != atoms_in_chunks {
        findings.push(
            Diagnostic::new(Code::E3401)
                .with_message("the residue table and the atom chunks disagree on the atom count")
                .with_context("residues account for", expected.to_string())
                .with_context("chunks hold", atoms_in_chunks.to_string()),
        );
    }
}

/// Every chain must reference an entity that exists.
fn check_entity_references(data: &StructureData, findings: &mut Diagnostics) {
    let topology = &data.topology;
    for chain in topology.chains.iter() {
        let Some(entity) = topology.chains.entity(chain) else {
            continue;
        };
        if entity.as_usize() >= topology.entities.len() {
            findings.push(
                Diagnostic::new(Code::E3005)
                    .with_context("chain", chain.to_string())
                    .with_context("entity", entity.to_string()),
            );
        }
    }
}

/// Models that share a topology must share an atom count.
fn check_coordinate_counts(data: &StructureData, findings: &mut Diagnostics) {
    let atoms = data.chunks.last().map_or(0, |chunk| chunk.atoms().end);
    let CoordinateStore::Dense { frames } = &data.coords else {
        if let CoordinateStore::Single(block) = &data.coords
            && block.len() != atoms
        {
            findings.push(
                Diagnostic::new(Code::E3401)
                    .with_context("positions", block.len().to_string())
                    .with_context("atoms", atoms.to_string()),
            );
        }
        return;
    };
    for (index, frame) in frames.iter().enumerate() {
        if frame.len() != atoms {
            findings.push(
                Diagnostic::new(Code::E3010)
                    .with_context("model", index.to_string())
                    .with_context("positions", frame.len().to_string())
                    .with_context("atoms", atoms.to_string()),
            );
        }
    }
}

/// Occupancy, where it was recorded, must lie between zero and one.
fn check_occupancies(data: &StructureData, findings: &mut Diagnostics) {
    for chunk in data.chunks.iter() {
        for local in 0..chunk.len() {
            let Some((occupancy, presence)) = chunk.occupancy(local) else {
                continue;
            };
            if presence.is_present() && !(0.0..=1.0).contains(&occupancy) {
                findings.push(
                    Diagnostic::new(Code::E3009)
                        .with_context("atom", (chunk.atoms().start + local).to_string())
                        .with_context("occupancy", occupancy.to_string()),
                );
            }
        }
    }
}

#[cfg(test)]
#[path = "validate_tests.rs"]
mod tests;
