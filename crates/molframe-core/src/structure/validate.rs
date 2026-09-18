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
use crate::Structure;
use crate::coords::CoordinateBlock;
use crate::diagnostic::{Code, Diagnostic, Diagnostics};
use crate::index::ResidueIndex;

/// Checks every invariant, returning what was violated.
///
/// An empty result means the structure is sound.
#[must_use]
pub fn validate(data: &StructureData) -> Vec<Diagnostic> {
    let mut findings = Diagnostics::new();
    let atoms = atoms_in_chunks(data);

    check_child_ranges(data, &mut findings);
    check_atom_coverage(data, atoms, &mut findings);
    check_entity_references(data, &mut findings);
    check_bonds(data, atoms, &mut findings);
    check_annotations(data, &mut findings);
    check_coordinate_counts(data, atoms, &mut findings);
    check_coordinates(data, &mut findings);
    check_occupancies(data, &mut findings);

    findings.finish()
}

/// Checks the invariants that a coordinate-only edit can change.
pub(crate) fn validate_coordinate_edit(data: &StructureData) -> Vec<Diagnostic> {
    let mut findings = Diagnostics::new();
    let atoms = atoms_in_chunks(data);
    check_coordinate_counts(data, atoms, &mut findings);
    check_coordinates(data, &mut findings);
    findings.finish()
}

/// Returns the atom extent declared by the structure's chunk sequence.
///
/// An empty chunk collection represents an empty atom extent.
fn atoms_in_chunks(data: &StructureData) -> u32 {
    data.chunks.last().map_or(0, |chunk| chunk.atoms().end)
}

/// Verifies that all annotation columns remain aligned with the atom table.
///
/// Every annotation must contain exactly one row for each atom.
fn check_annotations(data: &StructureData, findings: &mut Diagnostics) {
    let atoms = data.atom_count();

    for (name, column) in data.annotations.iter() {
        if column.len() == atoms {
            continue;
        }

        findings.push(
            Diagnostic::new(Code::E3011)
                .with_context("annotation", name)
                .with_context("column rows", column.len().to_string())
                .with_context("atoms", atoms.to_string()),
        );
    }
}

/// Verifies that every bond has two distinct, valid atom endpoints.
///
/// Self-bonds and endpoints beyond the atom extent are reported independently
/// so malformed bonds do not conceal additional structural violations.
fn check_bonds(data: &StructureData, atom_count: u32, findings: &mut Diagnostics) {
    for (position, bond) in data.bonds.iter().enumerate() {
        if bond.atom_a == bond.atom_b {
            findings.push(Diagnostic::new(Code::E3007).with_context("bond", position.to_string()));
        }

        check_bond_endpoint(position, bond.atom_a.get(), atom_count, findings);
        check_bond_endpoint(position, bond.atom_b.get(), atom_count, findings);
    }
}

/// Checks one bond endpoint against the current atom extent.
///
/// Invalid endpoints are reported with both the bond position and atom index.
fn check_bond_endpoint(bond: usize, atom: u32, atom_count: u32, findings: &mut Diagnostics) {
    if atom < atom_count {
        return;
    }

    findings.push(
        Diagnostic::new(Code::E3006)
            .with_context("bond", bond.to_string())
            .with_context("atom", atom.to_string()),
    );
}

/// Verifies coordinate values for every position declared present.
///
/// Ragged ensembles validate each independent model recursively. Single and
/// dense stores are traversed directly without repeated model lookups.
fn check_coordinates(data: &StructureData, findings: &mut Diagnostics) {
    match &data.coords {
        CoordinateStore::Single(block) => {
            check_coordinate_block(data, block, 0, findings);
        }
        CoordinateStore::Dense { frames } => {
            for (model, block) in frames.iter().enumerate() {
                let Ok(model) = u32::try_from(model) else {
                    findings
                        .push(Diagnostic::new(Code::E3001).with_context("model", "exceeds u32"));
                    return;
                };

                check_coordinate_block(data, block, model, findings);
            }
        }
        CoordinateStore::Ragged { models } => {
            check_ragged_models(models, findings);
        }
    }
}

/// Verifies finite coordinates for all atoms marked as positioned in one model.
///
/// Missing coordinates remain valid. Present coordinates must contain only
/// finite floating-point components.
fn check_coordinate_block(
    data: &StructureData,
    block: &CoordinateBlock,
    model: u32,
    findings: &mut Diagnostics,
) {
    let positions = block.as_slice();

    for chunk in data.chunks.iter() {
        let start = chunk.atoms().start;

        for local in 0..chunk.len() {
            if !chunk.has_position(local) {
                continue;
            }

            let Some(atom) = start.checked_add(local) else {
                findings.push(
                    Diagnostic::new(Code::E3001)
                        .with_context("model", model.to_string())
                        .with_context("atom", "index overflow"),
                );
                continue;
            };

            let Some(position) = positions.get(atom as usize) else {
                // The coordinate-count invariant reports truncated blocks.
                continue;
            };

            if position.iter().all(|value| value.is_finite()) {
                continue;
            }

            findings.push(
                Diagnostic::new(Code::E3008)
                    .with_context("model", model.to_string())
                    .with_context("atom", atom.to_string()),
            );
        }
    }
}

/// Recursively validates every independent structure in a ragged ensemble.
///
/// Diagnostics from each child carry the ragged-model index that produced
/// them, while nested ragged structures naturally accumulate model context.
fn check_ragged_models(models: &[Structure], findings: &mut Diagnostics) {
    for (model, structure) in models.iter().enumerate() {
        for finding in validate(structure.data()) {
            findings.push(finding.with_context("ragged model", model.to_string()));
        }
    }
}

/// Every parent's children must lie inside the child table, and siblings must
/// not overlap or run backwards.
fn check_child_ranges(data: &StructureData, findings: &mut Diagnostics) {
    let topology = &data.topology;
    let shared_topology = matches!(&data.coords, CoordinateStore::Dense { .. });

    let mut expected_chain = 0u32;
    let mut shared_range: Option<(u32, u32)> = None;

    for model in topology.models.iter() {
        let Some(range) = topology.models.chains(model) else {
            findings.push(
                Diagnostic::new(Code::E3002)
                    .with_message("a model has no chain range")
                    .with_context("model", model.to_string()),
            );
            continue;
        };

        if range.start > range.end {
            findings.push(
                Diagnostic::new(Code::E3003)
                    .with_message("a model's chain range runs backwards")
                    .with_context("model", model.to_string()),
            );
        }

        if range_exceeds_table(range.start, range.end, topology.chains.len()) {
            findings.push(
                Diagnostic::new(Code::E3002)
                    .with_message("a model's chains run past the chain table")
                    .with_context("model", model.to_string()),
            );
        }

        if shared_topology {
            let current = (range.start, range.end);

            match shared_range {
                Some(expected) if expected != current => {
                    findings.push(
                        Diagnostic::new(Code::E3003)
                            .with_message("dense models do not share one chain range")
                            .with_context("model", model.to_string()),
                    );
                }
                None => shared_range = Some(current),
                Some(_) => {}
            }
        } else {
            if range.start < expected_chain {
                findings.push(
                    Diagnostic::new(Code::E3003)
                        .with_message("models share chains")
                        .with_context("model", model.to_string()),
                );
            }

            expected_chain = expected_chain.max(range.end);
        }
    }

    let mut expected_residue = 0u32;

    for chain in topology.chains.iter() {
        let Some(range) = topology.chains.residues(chain) else {
            findings.push(
                Diagnostic::new(Code::E3002)
                    .with_message("a chain has no residue range")
                    .with_context("chain", chain.to_string()),
            );
            continue;
        };

        if range.start > range.end {
            findings.push(
                Diagnostic::new(Code::E3003)
                    .with_message("a chain's residue range runs backwards")
                    .with_context("chain", chain.to_string()),
            );
        }

        if range_exceeds_table(range.start, range.end, topology.residues.len()) {
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

        expected_residue = expected_residue.max(range.end);
    }
}

/// Determines whether either endpoint of a child range exceeds its table.
///
/// Both endpoints are checked so an inverted range whose `start` lies beyond
/// the table cannot evade validation merely because its `end` is smaller.
fn range_exceeds_table(start: u32, end: u32, table_len: usize) -> bool {
    let start = usize::try_from(start);
    let end = usize::try_from(end);

    match (start, end) {
        (Ok(start), Ok(end)) => start > table_len || end > table_len,
        _ => true,
    }
}

/// Every atom must belong to exactly one residue, which means the residues'
/// atom ranges must tile the atom order without gaps or overlaps.
fn check_atom_coverage(data: &StructureData, atoms_in_chunks: u32, findings: &mut Diagnostics) {
    let residues = &data.topology.residues;

    let Ok(residue_count) = u32::try_from(residues.len()) else {
        findings.push(Diagnostic::new(Code::E3001).with_context("residues", "exceeds u32"));
        return;
    };

    let mut expected = 0u32;

    for position in 0..residue_count {
        let residue = ResidueIndex::new(position);

        let Some(range) = residues.atoms(residue) else {
            findings.push(
                Diagnostic::new(Code::E3004)
                    .with_message("a residue has no atom range")
                    .with_context("residue", position.to_string()),
            );
            continue;
        };

        if range.start > range.end {
            findings.push(
                Diagnostic::new(Code::E3004)
                    .with_message("a residue's atom range runs backwards")
                    .with_context("residue", position.to_string()),
            );
        }

        if range.start != expected {
            findings.push(
                Diagnostic::new(Code::E3004)
                    .with_message("residues do not cover the atoms contiguously")
                    .with_context("residue", position.to_string())
                    .with_context("expected first atom", expected.to_string()),
            );
        }

        expected = expected.max(range.end);
    }

    if expected != atoms_in_chunks {
        findings.push(
            Diagnostic::new(Code::E3401)
                .with_message("the residue table and the atom chunks disagree on the atom count")
                .with_context("residues account for", expected.to_string())
                .with_context("chunks hold", atoms_in_chunks.to_string()),
        );
    }
}

/// Every chain must reference an entity that exists.
///
/// Missing entity references and references beyond the entity table are both
/// structural violations.
fn check_entity_references(data: &StructureData, findings: &mut Diagnostics) {
    let topology = &data.topology;

    for chain in topology.chains.iter() {
        let Some(entity) = topology.chains.entity(chain) else {
            findings.push(
                Diagnostic::new(Code::E3005)
                    .with_context("chain", chain.to_string())
                    .with_context("entity", "missing"),
            );
            continue;
        };

        if entity.as_usize() < topology.entities.len() {
            continue;
        }

        findings.push(
            Diagnostic::new(Code::E3005)
                .with_context("chain", chain.to_string())
                .with_context("entity", entity.to_string()),
        );
    }
}

/// Models that share a topology must share an atom count.
fn check_coordinate_counts(data: &StructureData, atoms: u32, findings: &mut Diagnostics) {
    match &data.coords {
        CoordinateStore::Single(block) => {
            if block.len() != atoms {
                findings.push(
                    Diagnostic::new(Code::E3401)
                        .with_context("positions", block.len().to_string())
                        .with_context("atoms", atoms.to_string()),
                );
            }
        }
        CoordinateStore::Dense { frames } => {
            for (index, frame) in frames.iter().enumerate() {
                if frame.len() == atoms {
                    continue;
                }

                findings.push(
                    Diagnostic::new(Code::E3010)
                        .with_context("model", index.to_string())
                        .with_context("positions", frame.len().to_string())
                        .with_context("atoms", atoms.to_string()),
                );
            }
        }
        CoordinateStore::Ragged { .. } => {}
    }
}

/// Occupancy, where it was recorded, must lie between zero and one.
fn check_occupancies(data: &StructureData, findings: &mut Diagnostics) {
    for chunk in data.chunks.iter() {
        let start = chunk.atoms().start;

        for local in 0..chunk.len() {
            let Some((occupancy, presence)) = chunk.occupancy(local) else {
                continue;
            };

            if !presence.is_present() || (0.0..=1.0).contains(&occupancy) {
                continue;
            }

            let Some(atom) = start.checked_add(local) else {
                findings.push(Diagnostic::new(Code::E3001).with_context("atom", "index overflow"));
                continue;
            };

            findings.push(
                Diagnostic::new(Code::E3009)
                    .with_context("atom", atom.to_string())
                    .with_context("occupancy", occupancy.to_string()),
            );
        }
    }
}

#[cfg(test)]
#[path = "validate_tests.rs"]
mod tests;
