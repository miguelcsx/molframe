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
//! The bonding reach comes from the process-wide `TableStore`: memory the
//! caller registered, else the profile and cache on disk, else the table the
//! `embedded-tables` feature compiles in. It is resolved once per read, not
//! once per model, and the table itself is decoded once per process.
//!
//! A ragged ensemble has no topology of its own: its atoms live in one
//! independent structure per model. Perception therefore runs on each model
//! and the shell is rebuilt around the enriched models, so a multi-model file
//! carries bonds and secondary structure exactly as a single-model one does.
//! A perception failure that the parse mode treats as an error ends the read
//! rather than being filed beside a structure that silently lacks its bonds.

use molframe_chem::{BondReach, Selected, TableError, TableStore};
use molframe_core::ExecutionContext;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics, Findings};
use molframe_core::io::ReadOptions;
use molframe_core::structure::{CoordinateStore, Structure as CoreStructure};

pub(super) fn enrich_read(
    result: Result<(CoreStructure, Vec<Diagnostic>), Findings>,
    options: &ReadOptions,
    context: &ExecutionContext,
) -> Result<(CoreStructure, Vec<Diagnostic>), Findings> {
    enrich_read_with(result, options, context, &TableStore::global())
}

/// [`enrich_read`] with the reach resolved from `store`.
pub(super) fn enrich_read_with(
    result: Result<(CoreStructure, Vec<Diagnostic>), Findings>,
    options: &ReadOptions,
    context: &ExecutionContext,
    store: &TableStore,
) -> Result<(CoreStructure, Vec<Diagnostic>), Findings> {
    let (structure, mut findings) = result?;
    if options.only_atomic_coords || options.perception == molframe_core::io::Perception::Skip {
        return Ok((structure, findings));
    }
    let reach = store.select_optional::<BondReach>("--bond-reach");
    let structure = match structure.ragged_models() {
        Some(models) => {
            let mut enriched = Vec::with_capacity(models.len());
            for model in models {
                enriched.push(enrich_model(model.clone(), &reach, &mut findings, context)?);
            }
            let mut data = structure.data().clone();
            data.coords = CoordinateStore::Ragged { models: enriched };
            CoreStructure::new(data)
        }
        None => enrich_model(structure, &reach, &mut findings, context)?,
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

/// The three stages on one model, with the bonding reach `reach` resolved to.
fn enrich_model(
    structure: CoreStructure,
    reach: &Result<Option<Selected<BondReach>>, TableError>,
    findings: &mut Vec<Diagnostic>,
    context: &ExecutionContext,
) -> Result<CoreStructure, Findings> {
    let structure = match reach {
        Ok(Some(reach)) => {
            match molframe_chem::perceive_bonds_in(&structure, &reach.provided.table, context) {
                Ok(structure) => structure,
                // A cancellation or an exhausted budget is the caller's own limit
                // taking effect, not a defect of the file: the read stops.
                Err(finding)
                    if matches!(finding.code(), Code::E1902 | Code::E1904 | Code::E7001) =>
                {
                    return Err(Findings::from(finding));
                }
                Err(finding) => {
                    findings.push(finding);
                    structure
                }
            }
        }
        Ok(None) => {
            findings.push(Diagnostic::new(Code::W3207));
            structure
        }
        Err(error) => {
            findings.push(Diagnostic::from(error.clone()));
            structure
        }
    };
    let structure = match molframe_chem::annotate_standard_components(&structure) {
        Ok(structure) => structure,
        Err(finding) => {
            findings.push(finding);
            structure
        }
    };
    let (states, sources) = merge_secondary(
        structure.data(),
        molframe_chem::assign_secondary_structure_in(&structure, context),
    );
    let mut data = structure.data().clone();
    data.secondary_structure = states.into();
    data.secondary_source = sources.into();
    Ok(CoreStructure::new(data))
}

/// Keeps, per residue, whichever assignment has the higher source rank, so a
/// deposited state always survives and analysis fills only what the file left
/// unassigned.
fn merge_secondary(
    data: &molframe_core::StructureData,
    inferred: Vec<molframe_core::SecondaryAssignment>,
) -> (
    Vec<molframe_core::SecondaryStructure>,
    Vec<molframe_core::SecondarySource>,
) {
    inferred
        .into_iter()
        .enumerate()
        .map(|(residue, inferred)| {
            let current = match (
                data.secondary_structure.get(residue),
                data.secondary_source.get(residue),
            ) {
                (Some(state), Some(source)) => molframe_core::SecondaryAssignment {
                    state: *state,
                    source: *source,
                },
                _ => molframe_core::SecondaryAssignment::default(),
            };
            let kept = if inferred.source.rank() > current.source.rank() {
                inferred
            } else {
                current
            };
            (kept.state, kept.source)
        })
        .unzip()
}
