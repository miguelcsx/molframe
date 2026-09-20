//! The rigid-transform verb over selected atoms.

use molframe_core::diagnostic::{Code, Diagnostic, Findings};
use molframe_core::selection::AtomSelection;
use molframe_geom::Rigid;

use crate::structure::Structure;

/// Applies one rigid transform to selected atoms in every dense model.
///
/// The coordinate transaction owns storage and generation tracking; the
/// geometric crate owns the transform arithmetic. This facade only joins the
/// two capabilities and publishes the resulting immutable snapshot.
///
/// # Errors
///
/// Returns a diagnostic when a selected atom is outside the topology, the
/// structure is a ragged ensemble, or the transformed snapshot is invalid.
pub fn transform(
    structure: &Structure,
    selection: &AtomSelection,
    rigid: &Rigid,
) -> Result<Structure, Findings> {
    let core = structure.engine();
    if core.ragged_models().is_some() {
        return Err(Diagnostic::new(Code::E6008).into());
    }
    if let Some(atom) = selection.iter().find(|atom| *atom >= core.atom_count()) {
        return Err(Diagnostic::new(Code::E6009)
            .with_context("atom", atom.to_string())
            .into());
    }

    let mut editor = core.edit();
    editor
        .transform(selection, |position| rigid.apply(position))
        .map_err(Findings::from)?;
    editor.commit().map(Structure::from).map_err(Findings::from)
}
