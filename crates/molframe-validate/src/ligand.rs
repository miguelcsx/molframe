//! Ligand bond geometry.
//!
//! The same covalent-radius sanity check the whole bond table gets, narrowed to
//! the bonds within het (non-polymer) residues — ligands, cofactors, modified
//! groups — where a wrong bond length most often hides. Standard-residue bonds
//! are left to the general check; here a bond is reported only when both its
//! atoms belong to het residues.

use molframe_core::index::AtomIndex;
use molframe_core::structure::{AtomRef, Structure};

use crate::covalent::covalent_deviation;
use crate::covalent_length::BondDeviation;

/// Ligand-bond outliers and explicit assessment coverage.
#[derive(Clone, Debug, PartialEq)]
pub struct LigandGeometryReport {
    /// Assessed bonds outside the requested tolerance.
    pub outliers: Vec<BondDeviation>,
    /// Bonds within hetero components selected for assessment.
    pub intended: usize,
    /// Selected bonds with elements and finite coordinates available.
    pub assessed: usize,
}

/// Assesses ligand bond geometry with explicit coverage.
#[must_use]
pub fn ligand_geometry(structure: &Structure, tolerance: f32) -> LigandGeometryReport {
    if !structure.data().bonds.is_available() {
        return LigandGeometryReport {
            outliers: Vec::new(),
            intended: 0,
            assessed: 0,
        };
    }
    let positions = structure.positions();
    let mut report = LigandGeometryReport {
        outliers: Vec::new(),
        intended: 0,
        assessed: 0,
    };
    for bond in structure.data().bonds.iter() {
        if !is_het(structure, bond.atom_a) || !is_het(structure, bond.atom_b) {
            continue;
        }
        report.intended += 1;
        let (Some(element_a), Some(element_b)) = (
            element_of(structure, bond.atom_a),
            element_of(structure, bond.atom_b),
        ) else {
            continue;
        };
        let (Some(&a), Some(&b)) = (
            positions.get(bond.atom_a.as_usize()),
            positions.get(bond.atom_b.as_usize()),
        ) else {
            continue;
        };
        let Some(deviation) = covalent_deviation(a, element_a, b, element_b) else {
            continue;
        };
        report.assessed += 1;
        if deviation.delta.abs() > tolerance {
            report.outliers.push(BondDeviation {
                atom_a: bond.atom_a,
                atom_b: bond.atom_b,
                observed: deviation.observed,
                expected: deviation.expected,
                deviation: deviation.delta,
            });
        }
    }
    report
}

/// Flags bonds within het residues whose length departs from the covalent-radii
/// sum by more than `tolerance`.
///
/// Returns an empty list when the structure has no bond table. Findings keep the
/// bond table's order.
///
/// Runs in `O(bonds)` time.
#[must_use]
pub fn ligand_geometry_outliers(structure: &Structure, tolerance: f32) -> Vec<BondDeviation> {
    ligand_geometry(structure, tolerance).outliers
}

/// Whether an atom belongs to a het (non-polymer) residue.
fn is_het(structure: &Structure, atom: AtomIndex) -> bool {
    match structure.data().atom(atom).and_then(AtomRef::residue) {
        Some(residue) => residue.is_het(),
        None => false,
    }
}

/// The element of an atom by index, if recorded.
fn element_of(structure: &Structure, atom: AtomIndex) -> Option<molframe_core::element::Element> {
    structure.data().atom(atom)?.element()
}

#[cfg(test)]
#[path = "ligand_tests.rs"]
mod tests;
