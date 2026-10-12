//! Per-residue solvent-accessible surface area and its largest possible value.
//!
//! The area is [`shrake_rupley`] summed over the atoms of each residue, so this
//! module only chooses the atoms and groups the result. The relative
//! accessibility a caller divides it by comes from [`max_residue_area`], the
//! Tien et al. (2013) theoretical maximum for each standard amino acid.
//!
//! What is in the calculation: the polymer residues (not `HETATM`) of the
//! structure's coordinates, hydrogens left out, one alternate conformer per
//! residue (the first label met). What is not: waters, ligands, ions and modified
//! residues written as `HETATM` neither carry an area nor shield their
//! neighbours, so a pocket a ligand fills reads as exposed. The counts of what
//! was left out are returned so a report can say so.

use molframe_chem::{RadiusSet, vdw_radius};
use molframe_core::ExecutionContext;
use molframe_core::structure::Structure;
use molframe_surface::{SasaError, shrake_rupley};

/// The table's name as a report records it.
pub const MAX_RESIDUE_AREA_REFERENCE: &str =
    "Tien et al. 2013, theoretical maximum ASA (Gly-X-Gly extended, all-atom)";

/// Where the table comes from.
pub const MAX_RESIDUE_AREA_CITATION: &str = "Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum allowed solvent accessibilities of residues in proteins. PLoS ONE 8(11): e80635";

/// The theoretical maximum accessible area of a standard amino acid, Å².
///
/// `None` for anything else (modified residues, unknown names), which has no
/// relative accessibility rather than a guessed one.
#[must_use]
pub fn max_residue_area(name: &str) -> Option<f64> {
    Some(match name {
        "ALA" => 129.0,
        "ARG" => 274.0,
        "ASN" => 195.0,
        "ASP" => 193.0,
        "CYS" => 167.0,
        "GLN" => 225.0,
        "GLU" => 223.0,
        "GLY" => 104.0,
        "ILE" => 197.0,
        "LEU" => 201.0,
        "LYS" => 236.0,
        "HIS" | "MET" => 224.0,
        "PHE" => 240.0,
        "PRO" => 159.0,
        "SER" => 155.0,
        "THR" => 172.0,
        "TRP" => 285.0,
        "TYR" => 263.0,
        "VAL" => 174.0,
        _ => return None,
    })
}

/// Why per-residue areas could not be computed.
#[derive(Debug, thiserror::Error)]
pub enum ResidueSasaError {
    /// No polymer atom with a known radius was left to measure.
    #[error("the structure has no polymer atom with a van der Waals radius to measure")]
    NoAtoms,
    /// The surface calculation was refused.
    #[error("the surface could not be computed: {0}")]
    Surface(#[from] SasaError),
}

/// The accessible area of one polymer residue.
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

/// The areas of every polymer residue, and what was left out.
#[derive(Clone, Debug, PartialEq)]
pub struct ResidueAreas {
    /// One row per residue, in structure order.
    pub residues: Vec<ResidueArea>,
    /// Hetero residues (waters, ligands, ions, modified residues) left out.
    pub excluded_hetero_residues: usize,
    /// Atoms of the kept residues left out: hydrogens, other conformers, no radius.
    pub skipped_atoms: usize,
}

/// The per-residue accessible area of a structure's polymer residues.
///
/// Radii are Bondi; hydrogens do not occlude. Runs one Shrake–Rupley pass over
/// the kept atoms, so it costs what [`shrake_rupley`] costs.
///
/// # Errors
///
/// [`ResidueSasaError::NoAtoms`] when nothing is left to measure;
/// [`ResidueSasaError::Surface`] when the surface calculation is refused (for
/// example, the memory budget).
pub fn residue_areas(
    structure: &Structure,
    probe: f32,
    points: u16,
    context: &ExecutionContext,
) -> Result<ResidueAreas, ResidueSasaError> {
    let positions = structure.positions();
    let mut coords = Vec::new();
    let mut radii = Vec::new();
    // (chain, number, name, index of its first kept atom, kept atom count)
    let mut spans: Vec<(String, i32, String, usize, usize)> = Vec::new();
    let (mut excluded, mut skipped) = (0, 0);
    for chain in structure.data().chains() {
        let label = chain
            .auth_label()
            .or_else(|| chain.label())
            .unwrap_or_default();
        for residue in chain.residues() {
            if residue.is_het() {
                excluded += 1;
                continue;
            }
            let number = residue
                .auth_seq_id()
                .or_else(|| residue.label_seq_id())
                .unwrap_or_default();
            let name = residue.name().unwrap_or("?");
            let start = coords.len();
            let mut conformer: Option<&str> = None;
            for atom in residue.atoms() {
                let radius = atom
                    .element()
                    .filter(|element| element.atomic_number() > 1)
                    .and_then(|element| vdw_radius(element, RadiusSet::Bondi));
                let (Some(radius), Some(at)) = (radius, positions.get(atom.index().get() as usize))
                else {
                    skipped += 1;
                    continue;
                };
                if let Some(alt) = atom.alt_label()
                    && *conformer.get_or_insert(alt) != alt
                {
                    skipped += 1;
                    continue;
                }
                coords.push(*at);
                radii.push(radius);
            }
            spans.push((
                label.to_owned(),
                number,
                name.to_owned(),
                start,
                coords.len() - start,
            ));
        }
    }
    if coords.is_empty() {
        return Err(ResidueSasaError::NoAtoms);
    }
    let areas = shrake_rupley(&coords, &radii, probe, points, context)?;
    let residues = spans
        .into_iter()
        .map(|(chain, number, name, start, count)| ResidueArea {
            chain,
            number,
            name,
            area: areas[start..start + count].iter().sum(),
            atoms: count,
        })
        .collect();
    Ok(ResidueAreas {
        residues,
        excluded_hetero_residues: excluded,
        skipped_atoms: skipped,
    })
}

#[cfg(test)]
#[path = "residue_sasa_tests.rs"]
mod tests;
