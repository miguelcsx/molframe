//! The perception pass that follows every file read.
//!
//! The stages run in this order on each model, each over the previous result:
//!
//! 1. bond perception (`molframe_chem::perceive_bonds_in`): file bonds, polymer
//!    links, then inferred covalent bonds from a parallel cell-grid search;
//!    this is the dominant stage, roughly a quarter of the parse time on a
//!    2.4-million-atom file;
//! 2. standard-component annotation: a table lookup per residue, negligible;
//! 3. secondary structure (DSSP): fills only residues the file left unassigned.
//!
//! `ReadOptions::perceive(false)` skips all three, as does `only_atomic_coords`
//! (which also drops the CIF extension projections). Nothing here is lazy: a
//! snapshot is immutable, so a stage is either part of the read or left out.
//!
//! A ragged ensemble has no topology of its own: its atoms live in one
//! independent structure per model. Perception therefore runs on each model
//! and the shell is rebuilt around the enriched models, so a multi-model file
//! carries bonds and secondary structure exactly as a single-model one does.
//! A perception failure that the parse mode treats as an error ends the read
//! rather than being filed beside a structure that silently lacks its bonds.

use super::enrich_model;
use molframe_core::ExecutionContext;
use molframe_core::diagnostic::{Diagnostic, Diagnostics, Findings};
use molframe_core::io::ReadOptions;
use molframe_core::structure::{CoordinateStore, Structure as CoreStructure};

pub(super) fn enrich_read(
    result: Result<(CoreStructure, Vec<Diagnostic>), Findings>,
    options: &ReadOptions,
    context: &ExecutionContext,
) -> Result<(CoreStructure, Vec<Diagnostic>), Findings> {
    let (structure, mut findings) = result?;
    if options.only_atomic_coords || options.perception == molframe_core::io::Perception::Skip {
        return Ok((structure, findings));
    }
    let structure = match structure.ragged_models() {
        Some(models) => {
            let mut enriched = Vec::with_capacity(models.len());
            for model in models {
                enriched.push(enrich_model(model.clone(), &mut findings, context)?);
            }
            let mut data = structure.data().clone();
            data.coords = CoordinateStore::Ragged { models: enriched };
            CoreStructure::new(data)
        }
        None => enrich_model(structure, &mut findings, context)?,
    };
    let mut ordered = Diagnostics::new();
    ordered.extend(findings);
    let failed = ordered.has_error(options.mode.strictness());
    let findings = ordered.finish();
    if failed {
        return Err(Findings::from(findings));
    }
    Ok((structure, findings))
}
