//! Per-residue solvent-accessible surface area and its largest possible value.
//!
//! The area is the per-atom area of [`molframe_surface::measure_atoms`] summed
//! over the atoms of each residue, so this module only chooses the atoms and
//! groups the result. Which residues take part, which radii they get and which
//! sampling method runs all come from one [`SasaOptions`], the same one the
//! facade, the command line and Python read. The relative accessibility a
//! caller divides the area by comes from [`max_residue_area`], the Tien et al.
//! (2013) theoretical maximum for each standard amino acid, or from a caller's own
//! [`MaxAreaTable`].
//!
//! What is in the calculation: the residues the options include (the polymer by
//! default, modified residues written as `HETATM` among it), hydrogens left
//! out, one alternate conformer per residue (the first label met). What is not:
//! excluded waters, ligands and ions neither carry an area nor shield their
//! neighbours, so a pocket a ligand fills reads as exposed unless
//! [`SasaOptions::hetero`] brings the ligand in. The counts of what was left
//! out are returned so a report can say so.

pub use crate::max_area_defaults::{MAX_RESIDUE_AREA_CITATION, MAX_RESIDUE_AREA_REFERENCE};
use crate::max_area_table::MaxAreaTable;
use molframe_chem::{AppliedRadiusPolicy, included_atoms};
use molframe_core::ExecutionContext;
use molframe_core::structure::Structure;
use molframe_surface::{SasaError, SasaOptions, measure_atoms};

/// The theoretical maximum accessible area of a standard amino acid, Å², from the built-in
/// table.
///
/// Protonation variants (`HID`, `CYX`, `ASH`, `GLH`, `LYN`, ...) take their
/// parent's value. `None` for anything else (modified residues, unknown names), which has no
/// relative accessibility rather than a guessed one. Pass a [`MaxAreaTable`] to use other values.
#[must_use]
pub fn max_residue_area(name: &str) -> Option<f64> {
    MaxAreaTable::default().max_area(name)
}

/// Why per-residue areas could not be computed.
#[derive(Debug, thiserror::Error)]
pub enum ResidueSasaError {
    /// No included atom with a known radius was left to measure.
    #[error("the structure has no included atom with a van der Waals radius to measure")]
    NoAtoms,
    /// The surface calculation was refused.
    #[error("the surface could not be computed: {0}")]
    Surface(#[from] SasaError),
}

/// The accessible area of one residue.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidueArea {
    /// Author chain label, else the label chain, else empty.
    pub chain: String,
    /// Author sequence number, else the label one, else zero.
    pub number: i32,
    /// Residue name, `?` when it has none.
    pub name: String,
    /// Solvent-accessible area, Å².
    pub area: f64,
    /// Heavy atoms that carried an area.
    pub atoms: usize,
}

/// The areas of every included residue, and what was left out.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidueAreas {
    /// One row per residue, in structure order.
    pub residues: Vec<ResidueArea>,
    /// Residues the heteroatom options left out.
    pub excluded_hetero_residues: usize,
    /// Atoms of the kept residues left out: hydrogens, other conformers, no radius.
    pub skipped_atoms: usize,
    /// The radius tables and choices the areas were computed with.
    pub policy: AppliedRadiusPolicy,
}

/// The per-residue accessible area of the residues `options` include.
///
/// Hydrogens do not occlude. Runs one pass over the kept atoms with the
/// method `options` names, so it costs what that kernel costs.
///
/// # Errors
///
/// [`ResidueSasaError::NoAtoms`] when nothing is left to measure;
/// [`ResidueSasaError::Surface`] when the surface calculation is refused (for
/// example, the memory budget).
pub fn residue_areas(
    structure: &Structure,
    options: &SasaOptions,
    context: &ExecutionContext,
) -> Result<ResidueAreas, ResidueSasaError> {
    let mut radii = options.atom_radii(structure);
    let included = included_atoms(structure, options.hetero);
    // (chain, number, name, atom indices that take part)
    let mut spans: Vec<(String, i32, String, Vec<usize>)> = Vec::new();
    let (mut excluded, mut skipped) = (0, 0);
    for chain in structure.data().chains() {
        let label = match chain.auth_label().or_else(|| chain.label()) {
            Some(label) => label,
            None => Default::default(),
        };
        for residue in chain.residues() {
            let first = residue.atoms().next().map(|atom| atom.index().as_usize());
            if !first.is_some_and(|atom| matches!(included.get(atom), Some(true))) {
                excluded += 1;
                continue;
            }
            let number = match residue.auth_seq_id().or_else(|| residue.label_seq_id()) {
                Some(number) => number,
                None => Default::default(),
            };
            let name = match residue.name() {
                Some(name) => name,
                None => "?",
            };
            let mut kept = Vec::new();
            let mut conformer: Option<&str> = None;
            for atom in residue.atoms() {
                let index = atom.index().as_usize();
                let hydrogen = atom
                    .element()
                    .is_some_and(|element| element.atomic_number() <= 1);
                let other_conformer = atom
                    .alt_label()
                    .is_some_and(|alt| *conformer.get_or_insert(alt) != alt);
                let radius = radii
                    .get(index)
                    .copied()
                    .filter(|radius| radius.is_finite());
                if hydrogen || other_conformer || radius.is_none() {
                    if let Some(slot) = radii.get_mut(index) {
                        *slot = f32::NAN;
                    }
                    skipped += 1;
                    continue;
                }
                kept.push(index);
            }
            spans.push((label.to_owned(), number, name.to_owned(), kept));
        }
    }
    if spans.iter().all(|(_, _, _, atoms)| atoms.is_empty()) {
        return Err(ResidueSasaError::NoAtoms);
    }
    let areas = measure_atoms(structure.positions(), &radii, options, context)?;
    let residues = spans
        .into_iter()
        .map(|(chain, number, name, atoms)| ResidueArea {
            chain,
            number,
            name,
            area: atoms.iter().map(|&atom| areas[atom]).sum(),
            atoms: atoms.len(),
        })
        .collect();
    Ok(ResidueAreas {
        residues,
        excluded_hetero_residues: excluded,
        skipped_atoms: skipped,
        policy: options.applied_policy(),
    })
}

#[cfg(test)]
#[path = "residue_sasa_tests.rs"]
mod tests;
