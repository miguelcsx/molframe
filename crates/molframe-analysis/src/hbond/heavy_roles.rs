//! Donor and acceptor typing from residue and atom names, for structures that
//! carry no hydrogens.
//!
//! Without a hydrogen a donor cannot be read from connectivity, so the role of
//! each heavy atom comes from what its residue is, following the conventions of
//! `McDonald` and Thornton's HBPLUS and of Baker and Hubbard:
//!
//! - protein backbone `N` donates (except proline, whose amide has no hydrogen,
//!   unless it is the chain's first residue), `O` and `OXT` accept;
//! - side chains: Arg `NE NH1 NH2`, Asn `ND2`, Gln `NE2`, Lys `NZ` and Trp `NE1`
//!   donate; Asn `OD1`, Gln `OE1`, Asp `OD1 OD2` and Glu `OE1 OE2` accept; His
//!   `ND1 NE2`, Ser `OG`, Thr `OG1` and Tyr `OH` do both;
//! - protonation variants follow their declared state: `HID`/`HSD` donate from
//!   `ND1` and accept at `NE2`, `HIE`/`HSE` the reverse, `HIP`/`HSP` donate from
//!   both; `ASH`/`GLH` accept at the carbonyl oxygen and do both at the hydroxyl;
//!   `LYN` `NZ` does both; `CYX`/`CYM` `SG` only accept;
//! - Cys `SG` (both) and Met `SD` (acceptor) are weak and only typed when asked;
//! - water oxygen does both;
//! - nucleic acids: base amino and imino nitrogens donate, ring nitrogens and
//!   carbonyl oxygens accept, phosphate and ester oxygens accept, `O2'` does both;
//! - every other group: a nitrogen or oxygen can do both, unless the structure
//!   carries donor and acceptor annotations, which are then authoritative.
//!
//! Typing costs one pass over the atoms.

use crate::protonation::parent_residue;
use molframe_core::column::Presence;
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use molframe_core::{AtomAnnotation, AtomIndex};

pub(super) const DONOR: u8 = 1;
pub(super) const ACCEPTOR: u8 = 2;
const BOTH: u8 = DONOR | ACCEPTOR;

/// The donor and acceptor atoms of a structure.
pub(super) struct Roles {
    pub(super) donors: AtomSelection,
    pub(super) acceptors: AtomSelection,
}

/// Types every atom. `include_sulfur` admits Cys `SG` and Met `SD`.
pub(super) fn assign_roles(structure: &Structure, include_sulfur: bool) -> Roles {
    let count = structure.atom_count();
    let mut bits = vec![0u8; count as usize];
    let annotated = annotated_roles(structure);
    for chain in structure.data().chains() {
        let mut first = true;
        for residue in chain.residues() {
            let starts_chain = std::mem::replace(&mut first, false);
            for atom in residue.atoms() {
                let (Some(component), Some(name)) = (atom.component_name(), atom.name()) else {
                    continue;
                };
                let index = atom.index().get() as usize;
                let role = standard_role(component, name, starts_chain, include_sulfur)
                    .or_else(|| group_role(structure, atom.index(), annotated));
                if let (Some(role), Some(slot)) = (role, bits.get_mut(index)) {
                    *slot = role;
                }
            }
        }
    }
    let select = |mask: u8| {
        AtomSelection::from_sorted(
            (0..count)
                .filter(|&atom| bits.get(atom as usize).is_some_and(|bits| bits & mask != 0))
                .collect(),
        )
    };
    Roles {
        donors: select(DONOR),
        acceptors: select(ACCEPTOR),
    }
}

/// Role of an atom in a standard polymer residue or water; `None` when the
/// residue is neither, so the caller falls back to the generic rule.
fn standard_role(component: &str, atom: &str, starts_chain: bool, sulfur: bool) -> Option<u8> {
    if is_water(component) {
        return Some(if atom.starts_with('O') { BOTH } else { 0 });
    }
    if is_amino_acid(component) {
        return Some(protein_role(component, atom, starts_chain, sulfur));
    }
    nucleic_role(component, atom)
}

fn is_water(component: &str) -> bool {
    matches!(component, "HOH" | "WAT" | "H2O" | "DOD")
}

fn is_amino_acid(component: &str) -> bool {
    matches!(
        parent_residue(component),
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
    )
}

fn protein_role(component: &str, atom: &str, starts_chain: bool, sulfur: bool) -> u8 {
    if let Some(role) = variant_role(component, atom, sulfur) {
        return role;
    }
    match (parent_residue(component), atom) {
        ("PRO", "N") => {
            if starts_chain {
                DONOR
            } else {
                0
            }
        }
        (_, "N")
        | ("ARG", "NE" | "NH1" | "NH2")
        | ("ASN", "ND2")
        | ("GLN", "NE2")
        | ("LYS", "NZ")
        | ("TRP", "NE1") => DONOR,
        (_, "O" | "OXT") => ACCEPTOR,
        ("ASN", "OD1") | ("GLN", "OE1") | ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2") => {
            ACCEPTOR
        }
        ("HIS", "ND1" | "NE2") | ("SER", "OG") | ("THR", "OG1") | ("TYR", "OH") => BOTH,
        ("CYS", "SG") if sulfur => BOTH,
        ("MET", "SD") if sulfur => ACCEPTOR,
        _ => 0,
    }
}

/// Side-chain atoms whose role depends on the protonation state the residue
/// name declares: a proton added makes an acceptor a donor, a proton removed
/// leaves only an acceptor. `None` for every other atom, which follows its
/// parent residue.
fn variant_role(component: &str, atom: &str, sulfur: bool) -> Option<u8> {
    Some(match (component, atom) {
        // A proton on the nitrogen makes it a donor, its absence an acceptor;
        // the doubly protonated ring donates from both.
        ("HID" | "HSD", "ND1") | ("HIE" | "HSE", "NE2") | ("HIP" | "HSP", "ND1" | "NE2") => DONOR,
        // Protonated carboxylic acid: the carbonyl oxygen accepts and the
        // hydroxyl oxygen does both. Neutral lysine amine: lone pair plus
        // hydrogens, so also both.
        ("ASH", "OD2") | ("GLH", "OE2") | ("LYN", "NZ") => BOTH,
        ("HID" | "HSD", "NE2") | ("HIE" | "HSE", "ND1") | ("ASH", "OD1") | ("GLH", "OE1") => {
            ACCEPTOR
        }
        // Disulphide and thiolate sulfur carry no hydrogen.
        ("CYX" | "CYM", "SG") => {
            if sulfur {
                ACCEPTOR
            } else {
                0
            }
        }
        _ => return None,
    })
}

fn nucleic_role(component: &str, atom: &str) -> Option<u8> {
    let base = match component {
        "DA" | "A" => Base::Adenine,
        "DC" | "C" => Base::Cytosine,
        "DG" | "G" => Base::Guanine,
        "DT" => Base::Thymine,
        "DU" | "U" => Base::Uracil,
        "DI" | "I" => Base::Hypoxanthine,
        _ => return None,
    };
    Some(match atom {
        "OP1" | "OP2" | "OP3" | "O1P" | "O2P" | "O3P" | "O5'" | "O3'" | "O4'" => ACCEPTOR,
        "O2'" => BOTH,
        _ => base.role(atom),
    })
}

#[derive(Clone, Copy)]
enum Base {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
    Uracil,
    Hypoxanthine,
}

impl Base {
    fn role(self, atom: &str) -> u8 {
        match (self, atom) {
            (Self::Adenine, "N6") | (Self::Cytosine, "N4") | (Self::Guanine, "N2") => DONOR,
            (Self::Guanine | Self::Hypoxanthine, "N1") | (Self::Thymine | Self::Uracil, "N3") => {
                DONOR
            }
            (Self::Adenine, "N1" | "N3" | "N7")
            | (Self::Cytosine, "N3" | "O2")
            | (Self::Guanine | Self::Hypoxanthine, "N3" | "N7" | "O6")
            | (Self::Thymine | Self::Uracil, "O2" | "O4") => ACCEPTOR,
            _ => 0,
        }
    }
}

type RoleColumns<'a> = (
    &'a molframe_core::AnnotationColumn<bool>,
    &'a molframe_core::AnnotationColumn<bool>,
);

fn annotated_roles(structure: &Structure) -> Option<RoleColumns<'_>> {
    let annotations = structure.annotations();
    let (Some(AtomAnnotation::Boolean(donors)), Some(AtomAnnotation::Boolean(acceptors))) = (
        annotations.get(molframe_core::HBOND_DONOR_ANNOTATION),
        annotations.get(molframe_core::HBOND_ACCEPTOR_ANNOTATION),
    ) else {
        return None;
    };
    let count = structure.atom_count();
    (donors.len() == count && acceptors.len() == count).then_some((donors, acceptors))
}

/// Generic rule for atoms outside the standard residues: the annotations when
/// present, otherwise any nitrogen or oxygen can do both.
fn group_role(
    structure: &Structure,
    atom: AtomIndex,
    annotated: Option<RoleColumns<'_>>,
) -> Option<u8> {
    if let Some((donors, acceptors)) = annotated {
        let flag = |column: &molframe_core::AnnotationColumn<bool>, bit: u8| {
            if matches!(column.get(atom.get()), Some((true, Presence::Present))) {
                bit
            } else {
                0
            }
        };
        return Some(flag(donors, DONOR) | flag(acceptors, ACCEPTOR));
    }
    let element = structure.data().atom(atom)?.element()?;
    Some(if matches!(element.atomic_number(), 7 | 8) {
        BOTH
    } else {
        0
    })
}
