//! Per-residue protein conformation: torsions, Ramachandran and rotamers.
//!
//! One pass over the residues, split into blocks that run in parallel and are
//! joined in block order, so the arrays are identical at any worker count. A
//! residue that cannot be evaluated (not a standard amino acid, missing atoms,
//! no neighbour across a chain break) carries `NaN` angles and the
//! [`NOT_EVALUATED`] class code instead of being dropped, which keeps every
//! array aligned with the structure's residue index.
//!
//! Complexity: O(residues) time, with constant work per residue (atom lookups
//! scan one residue's atoms; table lookups are a grid index). Parallel blocks of
//! 256 residues; the result columns are the only allocations, one `Vec` per
//! column plus one per block that is moved, not cloned, into the join.

use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

use super::BuiltinError;
use super::Layout;
use super::ReferenceSet;
use super::TableInfo;
use super::reference_data::RamaCategory;
use crate::numeric::f64_to_f32;

/// Class code for a residue the check does not apply to.
pub const NOT_EVALUATED: u8 = 255;
/// Residues per parallel block.
const BLOCK: usize = 256;

/// Peptide-bond planarity class, from the magnitude of omega.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OmegaClass {
    /// Within 30 degrees of planar trans.
    Trans,
    /// Within 30 degrees of planar cis.
    Cis,
    /// Neither: the peptide bond is twisted out of plane.
    Twisted,
}

impl OmegaClass {
    /// Classifies an omega angle in degrees.
    #[must_use]
    pub fn from_degrees(omega: f64) -> Self {
        let magnitude = omega.abs();
        if magnitude <= 30.0 {
            Self::Cis
        } else if magnitude >= 150.0 {
            Self::Trans
        } else {
            Self::Twisted
        }
    }

    /// Stable code used in the per-residue arrays.
    #[must_use]
    pub const fn code(self) -> u8 {
        match self {
            Self::Trans => 0,
            Self::Cis => 1,
            Self::Twisted => 2,
        }
    }
}

/// Per-residue protein conformation, one entry per structure residue.
#[derive(Clone, Debug, Default)]
pub struct ProteinTable {
    /// Phi in degrees, `NaN` when undefined.
    pub phi: Vec<f32>,
    /// Psi in degrees, `NaN` when undefined.
    pub psi: Vec<f32>,
    /// Omega of the bond to the previous residue, `NaN` when undefined.
    pub omega: Vec<f32>,
    /// Up to four chi angles per residue, row-major, `NaN` when undefined.
    pub chi: Vec<f32>,
    /// [`RamaCategory`] code of the distribution used.
    pub rama_category: Vec<u8>,
    /// Ramachandran density at the residue's (phi, psi).
    pub rama_density: Vec<f32>,
    /// 0 outlier, 1 allowed, 2 favoured, or [`NOT_EVALUATED`].
    pub rama_class: Vec<u8>,
    /// 0 outlier, 1 allowed, 2 favoured, or [`NOT_EVALUATED`].
    pub rotamer_class: Vec<u8>,
    /// [`OmegaClass`] code or [`NOT_EVALUATED`].
    pub omega_class: Vec<u8>,
    /// The reference tables the classes were scored against.
    pub tables: Vec<TableInfo>,
}

impl ProteinTable {
    /// Number of residues covered.
    #[must_use]
    pub fn len(&self) -> usize {
        self.phi.len()
    }

    /// Whether the table is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.phi.is_empty()
    }

    fn append(&mut self, other: Self) {
        self.phi.extend(other.phi);
        self.psi.extend(other.psi);
        self.omega.extend(other.omega);
        self.chi.extend(other.chi);
        self.rama_category.extend(other.rama_category);
        self.rama_density.extend(other.rama_density);
        self.rama_class.extend(other.rama_class);
        self.rotamer_class.extend(other.rotamer_class);
        self.omega_class.extend(other.omega_class);
    }

    fn push(&mut self, row: &Row) {
        self.phi.push(row.phi);
        self.psi.push(row.psi);
        self.omega.push(row.omega);
        self.chi.extend_from_slice(&row.chi);
        self.rama_category.push(row.rama_category);
        self.rama_density.push(row.rama_density);
        self.rama_class.push(row.rama_class);
        self.rotamer_class.push(row.rotamer_class);
        self.omega_class.push(row.omega_class);
    }
}

struct Row {
    phi: f32,
    psi: f32,
    omega: f32,
    chi: [f32; 4],
    rama_category: u8,
    rama_density: f32,
    rama_class: u8,
    rotamer_class: u8,
    omega_class: u8,
}

impl Row {
    const EMPTY: Self = Self {
        phi: f32::NAN,
        psi: f32::NAN,
        omega: f32::NAN,
        chi: [f32::NAN; 4],
        rama_category: NOT_EVALUATED,
        rama_density: f32::NAN,
        rama_class: NOT_EVALUATED,
        rotamer_class: NOT_EVALUATED,
        omega_class: NOT_EVALUATED,
    };
}

pub(super) fn torsion(layout: &Layout<'_>, atoms: [(usize, &str); 4]) -> Option<f64> {
    let [a, b, c, d] = atoms;
    molframe_geom::dihedral(
        layout.position(a.0, a.1)?,
        layout.position(b.0, b.1)?,
        layout.position(c.0, c.1)?,
        layout.position(d.0, d.1)?,
    )
    .map(f64::to_degrees)
}

pub(super) fn to_f32(value: Option<f64>) -> f32 {
    match value {
        Some(angle) => f64_to_f32(angle),
        None => f32::NAN,
    }
}

fn category(layout: &Layout<'_>, residue: usize, omega: Option<f64>) -> RamaCategory {
    let name = layout.name(residue);
    let before_proline = layout.name(residue + 1) == "PRO" && layout.peptide_linked(residue + 1);
    match name {
        "GLY" => RamaCategory::Gly,
        "PRO" if omega.is_some_and(|value| value.abs() < 90.0) => RamaCategory::CisPro,
        "PRO" => RamaCategory::TransPro,
        _ if before_proline => RamaCategory::PrePro,
        "ILE" | "VAL" => RamaCategory::IleVal,
        _ => RamaCategory::General,
    }
}

fn row(layout: &Layout<'_>, residue: usize, reference: &ReferenceSet) -> Row {
    let mut out = Row::EMPTY;
    if !layout.is_amino_acid(residue) {
        return out;
    }
    let previous = residue.wrapping_sub(1);
    let linked_before = layout.peptide_linked(residue);
    let linked_after = layout.peptide_linked(residue + 1);
    let omega = linked_before
        .then(|| {
            torsion(
                layout,
                [
                    (previous, "CA"),
                    (previous, "C"),
                    (residue, "N"),
                    (residue, "CA"),
                ],
            )
        })
        .flatten();
    let phi = linked_before
        .then(|| {
            torsion(
                layout,
                [
                    (previous, "C"),
                    (residue, "N"),
                    (residue, "CA"),
                    (residue, "C"),
                ],
            )
        })
        .flatten();
    let psi = linked_after
        .then(|| {
            torsion(
                layout,
                [
                    (residue, "N"),
                    (residue, "CA"),
                    (residue, "C"),
                    (residue + 1, "N"),
                ],
            )
        })
        .flatten();
    out.phi = to_f32(phi);
    out.psi = to_f32(psi);
    out.omega = to_f32(omega);
    if let Some(value) = omega {
        out.omega_class = OmegaClass::from_degrees(value).code();
    }
    if let (Some(phi), Some(psi)) = (phi, psi) {
        let kind = category(layout, residue, omega);
        let density = reference.rama.density(kind, phi, psi);
        out.rama_category = kind.code();
        out.rama_density = density;
        out.rama_class = reference.rama.classify(kind, density);
    }
    chi_and_rotamer(layout, residue, reference, &mut out);
    out
}

fn chi_and_rotamer(layout: &Layout<'_>, residue: usize, reference: &ReferenceSet, out: &mut Row) {
    let Some(definition) = reference.chi.get(layout.name(residue)) else {
        return;
    };
    let mut angles = [None; 4];
    let mut count = 0;
    for ((slot, names), angle) in out.chi.iter_mut().zip(&definition.chi).zip(&mut angles) {
        let value = torsion(
            layout,
            [
                (residue, names[0].as_str()),
                (residue, names[1].as_str()),
                (residue, names[2].as_str()),
                (residue, names[3].as_str()),
            ],
        );
        *slot = to_f32(value);
        *angle = value;
        count += 1;
    }
    let Some(table) = reference.rotamers.get(&definition.table) else {
        return;
    };
    let mut chi = [0.0f64; 4];
    for (slot, angle) in chi.iter_mut().zip(&angles[..count]) {
        let Some(angle) = angle else {
            return;
        };
        *slot = *angle;
    }
    if let Some(class) = table.classify(&chi[..count]) {
        out.rotamer_class = class.code();
    }
}

/// Computes the per-residue protein conformation of a laid-out structure.
///
/// Every threshold and table comes from `reference`; the result records which.
///
/// Runs in `O(residues)` time split across the context's workers.
///
/// # Errors
///
/// Returns [`BuiltinError`] if a worker panics.
pub fn protein_table(
    layout: &Layout<'_>,
    reference: &ReferenceSet,
    context: &ExecutionContext,
) -> Result<ProteinTable, BuiltinError> {
    let plan = BlockPlan::new(layout.len(), BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range| {
        let mut part = ProteinTable::default();
        for residue in range {
            part.push(&row(layout, residue, reference));
        }
        part
    })?;
    let mut table = ProteinTable {
        tables: reference.infos(),
        ..ProteinTable::default()
    };
    for block in blocks {
        table.append(block);
    }
    Ok(table)
}

#[cfg(test)]
#[path = "protein_tests.rs"]
mod tests;
