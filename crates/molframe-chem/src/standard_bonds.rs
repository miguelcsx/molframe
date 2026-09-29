//! Bond orders inside the standard polymer components.
//!
//! Geometry says two atoms are bonded; it cannot say how. For the standard
//! amino acids and nucleotides the order follows from IUPAC atom names alone,
//! so this table answers it without a component dictionary. Every lookup is a
//! match on two short names — `O(1)` time, no allocation — and a bond the table
//! does not name is single.
//!
//! Aromatic rings are reported as [`BondOrder::Aromatic`] rather than as a
//! Kekulé pattern: a renderer draws the delocalised ring, and a Kekulé choice
//! would be one arbitrary structure of several.

use molframe_core::BondOrder;

/// Order of the bond between `atom_a` and `atom_b` inside one residue of
/// `component`, for the standard polymer components; single otherwise.
#[must_use]
pub fn standard_bond_order(component: &str, atom_a: &str, atom_b: &str) -> BondOrder {
    let (low, high) = if atom_a <= atom_b {
        (atom_a, atom_b)
    } else {
        (atom_b, atom_a)
    };
    if in_ring(aromatic_ring(component), low) && in_ring(aromatic_ring(component), high) {
        return BondOrder::Aromatic;
    }
    if is_double(component, low, high) {
        BondOrder::Double
    } else {
        BondOrder::Single
    }
}

/// Whether `component` names water.
#[must_use]
pub fn is_water_component(component: &str) -> bool {
    matches!(
        component,
        "HOH" | "WAT" | "DOD" | "H2O" | "SOL" | "TIP3" | "TIP"
    )
}

/// Whether `component` is one of the standard or common amino acids whose
/// backbone carbonyl is named `C`/`O`.
#[must_use]
pub fn is_amino_acid_component(component: &str) -> bool {
    matches!(
        component,
        "ALA"
            | "ARG"
            | "ASN"
            | "ASP"
            | "CYS"
            | "GLN"
            | "GLU"
            | "GLY"
            | "HIS"
            | "ILE"
            | "LEU"
            | "LYS"
            | "MET"
            | "PHE"
            | "PRO"
            | "SER"
            | "THR"
            | "TRP"
            | "TYR"
            | "VAL"
            | "MSE"
            | "SEC"
            | "PYL"
            | "HID"
            | "HIE"
            | "HIP"
            | "HSD"
            | "HSE"
            | "HSP"
            | "UNK"
    )
}

/// Whether `component` is a standard ribo- or deoxyribonucleotide.
#[must_use]
pub fn is_nucleotide_component(component: &str) -> bool {
    matches!(
        component,
        "A" | "C" | "G" | "U" | "I" | "DA" | "DC" | "DG" | "DT" | "DU" | "DI"
    )
}

const PHENYL: &[&str] = &["CD1", "CD2", "CE1", "CE2", "CG", "CZ"];
const INDOLE: &[&str] = &["CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3", "NE1"];
const IMIDAZOLE: &[&str] = &["CD2", "CE1", "CG", "ND1", "NE2"];
const PURINE: &[&str] = &["C2", "C4", "C5", "C6", "C8", "N1", "N3", "N7", "N9"];
const PYRIMIDINE: &[&str] = &["C2", "C4", "C5", "C6", "N1", "N3"];

fn aromatic_ring(component: &str) -> &'static [&'static str] {
    match component {
        "PHE" | "TYR" => PHENYL,
        "TRP" => INDOLE,
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => IMIDAZOLE,
        "A" | "G" | "I" | "DA" | "DG" | "DI" => PURINE,
        "C" | "U" | "DC" | "DT" | "DU" => PYRIMIDINE,
        _ => &[],
    }
}

fn in_ring(ring: &[&str], atom: &str) -> bool {
    ring.contains(&atom)
}

/// Exocyclic and side-chain double bonds, with `low <= high`.
fn is_double(component: &str, low: &str, high: &str) -> bool {
    if is_amino_acid_component(component) && low == "C" && high == "O" {
        return true;
    }
    if is_nucleotide_component(component) && low == "OP1" && high == "P" {
        return true;
    }
    matches!(
        (component, low, high),
        ("ARG", "CZ", "NH2")
            | ("ASN" | "ASP", "CG", "OD1")
            | ("GLN" | "GLU", "CD", "OE1")
            | ("G" | "DG" | "I" | "DI", "C6", "O6")
            | ("C" | "DC" | "U" | "DU" | "DT", "C2", "O2")
            | ("U" | "DU" | "DT", "C4", "O4")
    )
}

#[cfg(test)]
#[path = "standard_bonds_tests.rs"]
mod tests;
