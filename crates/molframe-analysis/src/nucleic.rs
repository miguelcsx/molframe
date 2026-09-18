//! Nucleic backbone and glycosidic torsions from explicit semantic atom roles.

use molframe_chem::PolymerAtomRole;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::{AtomRef, ResidueRef, Structure};
use molframe_core::{AtomAnnotation, Presence};

/// The seven torsions of one nucleotide, in degrees where defined.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NucleicTorsions {
    /// The residue described.
    pub residue: ResidueIndex,
    /// α: O3(i−1)–P–O5–C5.
    pub alpha: Option<f64>,
    /// β: P–O5–C5–C4.
    pub beta: Option<f64>,
    /// γ: O5–C5–C4–C3.
    pub gamma: Option<f64>,
    /// δ: C5–C4–C3–O3.
    pub delta: Option<f64>,
    /// ε: C4–C3–O3–P(i+1).
    pub epsilon: Option<f64>,
    /// ζ: C3–O3–P(i+1)–O5(i+1).
    pub zeta: Option<f64>,
    /// χ: O4–C1–glycosidic-base–base-reference.
    pub chi: Option<f64>,
}

/// Why nucleic torsions could not be projected from a structure.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum NucleicTorsionError {
    /// An explicit typed polymer atom-role column is required.
    #[error("nucleic torsions require explicit CCD polymer atom-role annotations")]
    MissingRoles,
    /// More than one atom carries a role that must be unique in a residue.
    #[error("residue {residue} has multiple atoms for polymer role {role}")]
    AmbiguousRole {
        /// Affected residue.
        residue: ResidueIndex,
        /// Stable role bit code.
        role: i64,
    },
}

/// Computes backbone and glycosidic torsions from explicit atom roles and bonds.
///
/// Inter-residue continuity is an actual O3-to-phosphate bond. No atom-name
/// aliases, residue tables, or distance cutoff are embedded in this kernel.
///
/// # Errors
///
/// Returns [`NucleicTorsionError`] when semantic roles are absent or ambiguous.
pub fn nucleic_torsions(
    structure: &Structure,
) -> Result<Vec<NucleicTorsions>, NucleicTorsionError> {
    require_roles(structure)?;
    let mut records = Vec::new();
    for chain in structure.data().chains() {
        let residues: Vec<ResidueRef<'_>> = chain
            .residues()
            .filter(|residue| {
                residue_has_role(
                    structure,
                    *residue,
                    PolymerAtomRole::NUCLEIC_BACKBONE
                        .union(PolymerAtomRole::NUCLEIC_SUGAR)
                        .union(PolymerAtomRole::NUCLEIC_BASE_GROUP),
                )
            })
            .collect();
        for position in 0..residues.len() {
            let current = residues[position];
            let previous = position
                .checked_sub(1)
                .and_then(|index| residues.get(index).copied());
            let next = residues.get(position + 1).copied();
            records.push(torsions_of(structure, current, previous, next)?);
        }
    }
    records.sort_by_key(|record| record.residue.get());
    Ok(records)
}

fn torsions_of(
    structure: &Structure,
    current: ResidueRef<'_>,
    previous: Option<ResidueRef<'_>>,
    next: Option<ResidueRef<'_>>,
) -> Result<NucleicTorsions, NucleicTorsionError> {
    let previous = previous
        .map(|residue| linked(structure, residue, current).map(|linked| linked.then_some(residue)))
        .transpose()?
        .flatten();
    let next = next
        .map(|residue| linked(structure, current, residue).map(|linked| linked.then_some(residue)))
        .transpose()?
        .flatten();
    let current_atoms = NucleotideAtoms::project(structure, current)?;
    let previous_o3 = previous
        .map(|residue| role_position(structure, residue, PolymerAtomRole::NUCLEIC_O3))
        .transpose()?
        .flatten();
    let next_atoms = next
        .map(|residue| NucleotideAtoms::project(structure, residue))
        .transpose()?;
    Ok(NucleicTorsions {
        residue: current.index(),
        alpha: torsion(
            previous_o3,
            current_atoms.p,
            current_atoms.o5,
            current_atoms.c5,
        ),
        beta: torsion(
            current_atoms.p,
            current_atoms.o5,
            current_atoms.c5,
            current_atoms.c4,
        ),
        gamma: torsion(
            current_atoms.o5,
            current_atoms.c5,
            current_atoms.c4,
            current_atoms.c3,
        ),
        delta: torsion(
            current_atoms.c5,
            current_atoms.c4,
            current_atoms.c3,
            current_atoms.o3,
        ),
        epsilon: torsion(
            current_atoms.c4,
            current_atoms.c3,
            current_atoms.o3,
            next_atoms.and_then(|atoms| atoms.p),
        ),
        zeta: torsion(
            current_atoms.c3,
            current_atoms.o3,
            next_atoms.and_then(|atoms| atoms.p),
            next_atoms.and_then(|atoms| atoms.o5),
        ),
        chi: torsion(
            current_atoms.o4,
            current_atoms.c1,
            current_atoms.glycosidic,
            current_atoms.base_reference,
        ),
    })
}

#[derive(Clone, Copy)]
struct NucleotideAtoms {
    p: Option<[f32; 3]>,
    o5: Option<[f32; 3]>,
    c5: Option<[f32; 3]>,
    c4: Option<[f32; 3]>,
    c3: Option<[f32; 3]>,
    o3: Option<[f32; 3]>,
    o4: Option<[f32; 3]>,
    c1: Option<[f32; 3]>,
    glycosidic: Option<[f32; 3]>,
    base_reference: Option<[f32; 3]>,
}

impl NucleotideAtoms {
    fn project(
        structure: &Structure,
        residue: ResidueRef<'_>,
    ) -> Result<Self, NucleicTorsionError> {
        Ok(Self {
            p: role_position(structure, residue, PolymerAtomRole::NUCLEIC_PHOSPHATE)?,
            o5: role_position(structure, residue, PolymerAtomRole::NUCLEIC_O5)?,
            c5: role_position(structure, residue, PolymerAtomRole::NUCLEIC_C5)?,
            c4: role_position(structure, residue, PolymerAtomRole::NUCLEIC_C4)?,
            c3: role_position(structure, residue, PolymerAtomRole::NUCLEIC_C3)?,
            o3: role_position(structure, residue, PolymerAtomRole::NUCLEIC_O3)?,
            o4: role_position(structure, residue, PolymerAtomRole::NUCLEIC_O4)?,
            c1: role_position(structure, residue, PolymerAtomRole::NUCLEIC_C1)?,
            glycosidic: role_position(structure, residue, PolymerAtomRole::NUCLEIC_GLYCOSIDIC)?,
            base_reference: role_position(
                structure,
                residue,
                PolymerAtomRole::NUCLEIC_BASE_REFERENCE,
            )?,
        })
    }
}

fn linked(
    structure: &Structure,
    first: ResidueRef<'_>,
    second: ResidueRef<'_>,
) -> Result<bool, NucleicTorsionError> {
    let first_o3 = role_atom(structure, first, PolymerAtomRole::NUCLEIC_O3)?;
    let second_p = role_atom(structure, second, PolymerAtomRole::NUCLEIC_PHOSPHATE)?;
    let (Some(first_o3), Some(second_p)) = (first_o3, second_p) else {
        return Ok(false);
    };
    if !structure.data().bonds.is_available() {
        return Ok(false);
    }
    Ok(structure
        .data()
        .bonds
        .adjacency(structure.atom_count())
        .neighbours(first_o3.index())
        .binary_search(&second_p.index())
        .is_ok())
}

fn require_roles(structure: &Structure) -> Result<(), NucleicTorsionError> {
    if matches!(
        structure
            .annotations()
            .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION),
        Some(AtomAnnotation::Integer(_))
    ) {
        Ok(())
    } else {
        Err(NucleicTorsionError::MissingRoles)
    }
}

fn residue_has_role(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
) -> bool {
    residue
        .atoms()
        .any(|atom| atom_role(structure, atom).is_some_and(|role| role.intersects(required)))
}

fn role_position(
    structure: &Structure,
    residue: ResidueRef<'_>,
    required: PolymerAtomRole,
) -> Result<Option<[f32; 3]>, NucleicTorsionError> {
    Ok(role_atom(structure, residue, required)?.and_then(AtomRef::position))
}

fn role_atom<'a>(
    structure: &'a Structure,
    residue: ResidueRef<'a>,
    required: PolymerAtomRole,
) -> Result<Option<AtomRef<'a>>, NucleicTorsionError> {
    let mut matches = residue
        .atoms()
        .filter(|atom| atom_role(structure, *atom).is_some_and(|role| role.intersects(required)));
    let first = matches.next();
    if matches.next().is_some() {
        Err(NucleicTorsionError::AmbiguousRole {
            residue: residue.index(),
            role: required.code(),
        })
    } else {
        Ok(first)
    }
}

fn atom_role(structure: &Structure, atom: AtomRef<'_>) -> Option<PolymerAtomRole> {
    let AtomAnnotation::Integer(column) = structure
        .annotations()
        .get(molframe_core::POLYMER_ATOM_ROLE_ANNOTATION)?
    else {
        return None;
    };
    column
        .get(atom.index().get())
        .filter(|(_, presence)| *presence == Presence::Present)
        .and_then(|(code, _)| PolymerAtomRole::from_code(code))
}

fn torsion(
    first: Option<[f32; 3]>,
    second: Option<[f32; 3]>,
    third: Option<[f32; 3]>,
    fourth: Option<[f32; 3]>,
) -> Option<f64> {
    molframe_geom::dihedral(first?, second?, third?, fourth?).map(molframe_geom::degrees)
}

#[cfg(test)]
#[path = "nucleic_tests.rs"]
mod tests;
