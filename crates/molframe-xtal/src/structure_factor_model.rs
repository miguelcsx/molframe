//! Structure factors of a loaded structure from its own cell and symmetry.

use crate::{
    CellTransform, Complex64, Displacement, ScatteringSite, StructureFactorCalculator,
    StructureFactorError, SymmetryExt,
};
use molframe_core::index::{AtomIndex, ModelIndex};
use molframe_core::structure::Structure;
use molframe_core::{Code, Diagnostic};

/// Computes the X-ray structure factors of the first model at `reflections`.
///
/// Every atom row contributes with its occupancy and either its anisotropic
/// tensor or its isotropic `B`; the summation runs over the structure's full
/// space-group operation set, so the model is the asymmetric unit.
///
/// # Errors
///
/// Returns a diagnostic when the structure has no usable cell, no symmetry
/// operations, an atom without an element or form factor, or no first model.
pub fn structure_factors(
    structure: &Structure,
    reflections: &[[i32; 3]],
) -> Result<Vec<Complex64>, Diagnostic> {
    let cell = structure
        .data()
        .cell
        .filter(|cell| !cell.is_placeholder())
        .ok_or_else(|| Diagnostic::new(Code::E5004).with_context("cell", "absent"))?;
    let transform = CellTransform::new(&cell)?;
    let symmetry = structure
        .symmetry_set()
        .filter(|set| !set.operations().is_empty())
        .ok_or_else(|| Diagnostic::new(Code::E6016).with_context("symmetry", "no operations"))?;
    let positions = structure
        .model_positions(ModelIndex::new(0))
        .ok_or_else(|| Diagnostic::new(Code::E6003).with_context("model", "0"))?;
    let mut sites = Vec::with_capacity(positions.len());
    for (row, position) in positions.iter().enumerate() {
        let Ok(atom_number) = u32::try_from(row) else {
            break;
        };
        let index = AtomIndex::new(atom_number);
        let Some(atom) = structure.data().atom(index) else {
            continue;
        };
        let element = atom.element().ok_or_else(|| {
            Diagnostic::new(Code::E6016).with_context("atom", format!("{row} has no element"))
        })?;
        let displacement = match structure.data().anisotropy.for_atom(index) {
            Some(tensor) => Displacement::Anisotropic(tensor.map(f64::from)),
            // No recorded B means no displacement is modelled.
            None => Displacement::Isotropic(match atom.b_factor() {
                Some(b_factor) => f64::from(b_factor),
                None => 0.0,
            }),
        };
        sites.push(ScatteringSite {
            element,
            position: transform.to_fractional(position.map(f64::from)),
            // No recorded occupancy is the PDB convention for a full site.
            occupancy: match atom.occupancy() {
                Some(occupancy) => f64::from(occupancy),
                None => 1.0,
            },
            displacement,
        });
    }
    let calculator = StructureFactorCalculator::new(transform, symmetry.operations())
        .map_err(factor_diagnostic)?;
    calculator
        .calculate_many(&sites, reflections)
        .map_err(factor_diagnostic)
}

fn factor_diagnostic(error: StructureFactorError) -> Diagnostic {
    Diagnostic::new(Code::E6016).with_context("structure_factor", error.to_string())
}
