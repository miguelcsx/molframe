//! Guards for writing a multi-model structure.

use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::Select;
use molframe_core::structure::Structure;

/// Refuses a chain, residue or atom filter on a ragged ensemble that writes
/// more than one model, because its indices only describe the first model.
pub(super) fn refuse_foreign_filter(
    structure: &Structure,
    select: &impl Select,
) -> Result<(), Vec<Diagnostic>> {
    let Some(first) = structure.ragged_models().and_then(<[_]>::first) else {
        return Ok(());
    };
    let data = first.data();
    let narrows = data.chains().any(|chain| {
        !select.accept_chain(chain.index())
            || chain.residues().any(|residue| {
                !select.accept_residue(residue.index())
                    || residue
                        .atoms()
                        .any(|atom| !select.accept_atom(atom.index()))
            })
    });
    if narrows {
        return Err(vec![
            Diagnostic::new(Code::E6003)
                .with_context(
                    "reason",
                    "ragged ensemble: an index filter cannot span models with different topology; select one model",
                )
                .with_context("models", structure.model_count().to_string()),
        ]);
    }
    Ok(())
}

/// Refuses a later model whose coordinates the fixed columns cannot hold, which
/// the first-model capacity check cannot see.
pub(super) fn check_frame(
    structure: &Structure,
    frame: &[[f32; 3]],
    model: usize,
    select: &impl Select,
) -> Result<(), Vec<Diagnostic>> {
    let mut refusals = Vec::new();
    for atom in structure.data().atoms() {
        if !select.accept_atom(atom.index()) {
            continue;
        }
        let Some(position) = frame.get(atom.index().as_usize()) else {
            continue;
        };
        if position
            .iter()
            .any(|value| !super::capacity::coordinate_fits(f64::from(*value)))
        {
            refusals.push(
                Diagnostic::new(Code::E4104)
                    .with_context("model", (model + 1).to_string())
                    .with_context("atom", atom.index().to_string())
                    .with_context("position", format!("{position:?}")),
            );
        }
    }
    if refusals.is_empty() {
        Ok(())
    } else {
        Err(refusals)
    }
}
