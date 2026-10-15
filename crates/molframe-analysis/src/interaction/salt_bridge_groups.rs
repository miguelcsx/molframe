//! Charged groups for salt-bridge detection.
//!
//! The Chemical Component Dictionary records the neutral protonation state of
//! the standard amino acids, so a carboxylate oxygen or a lysine nitrogen
//! carries formal charge zero there even though both are charged at
//! physiological pH. Membership therefore comes from two sources, joined:
//!
//! 1. the formal-charge annotation, which is authoritative for ligands and any
//!    atom it marks as charged;
//! 2. the titratable groups of the standard amino acids: the Asp and Glu
//!    carboxylate oxygens, the Lys amine nitrogen, the Arg guanidinium
//!    nitrogens, and the chain termini (C-terminal `OXT`, and the backbone `N`
//!    of a chain's first residue). Histidine is not treated as charged: its
//!    protonation depends on the local environment, which a bare structure
//!    does not record.
//!
//! This is the side-chain model of Barlow and Thornton: nitrogen and oxygen
//! atoms of oppositely charged groups within 4 Å.

use molframe_core::structure::Structure;

/// The charge sign a titratable atom carries at physiological pH.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Sign {
    Anion,
    Cation,
}

const AMINO_ACIDS: [&str; 20] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

fn side_chain_sign(component: &str, atom: &str) -> Option<Sign> {
    match (component, atom) {
        ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2") => Some(Sign::Anion),
        ("LYS", "NZ") | ("ARG", "NE" | "NH1" | "NH2") => Some(Sign::Cation),
        _ => None,
    }
}

/// Every atom that is a titratable charged group, with its sign, in atom order.
pub(crate) fn titratable_atoms(structure: &Structure) -> Vec<(u32, Sign)> {
    let mut found = Vec::new();
    for chain in structure.data().chains() {
        let mut first = true;
        for residue in chain.residues() {
            let starts_chain = std::mem::replace(&mut first, false);
            for atom in residue.atoms() {
                let (Some(component), Some(name)) = (atom.component_name(), atom.name()) else {
                    continue;
                };
                let sign = if !AMINO_ACIDS.contains(&component) {
                    None
                } else if name == "OXT" {
                    Some(Sign::Anion)
                } else if starts_chain && name == "N" {
                    Some(Sign::Cation)
                } else {
                    side_chain_sign(component, name)
                };
                if let Some(sign) = sign {
                    found.push((atom.index().get(), sign));
                }
            }
        }
    }
    found.sort_unstable_by_key(|&(atom, _)| atom);
    found
}

/// Whether the structure has any atom in a titratable charged group.
pub(crate) fn has_titratable_groups(structure: &Structure) -> bool {
    !titratable_atoms(structure).is_empty()
}
