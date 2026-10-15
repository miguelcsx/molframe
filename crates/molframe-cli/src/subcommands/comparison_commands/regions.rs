//! RMSD over an explicitly selected region or ligand of corresponding atoms.
//!
//! The model shares the reference's atom numbering, so a selection made in the
//! reference names the same atoms in the model.

use super::command::ComparisonOptions;
use crate::MetricChoice;
use crate::args::AlignmentChoice;
use crate::exit::Exit;
use crate::report::Context;
use molframe::compare::workflow::{ComparisonAlignment, PointMapping, PointMatch, align_mapping};

type Rows = Result<Vec<(&'static str, f64)>, Exit>;

fn selected(
    reference: &molframe::Structure,
    query: &str,
    context: Context,
) -> Result<Vec<usize>, Exit> {
    let evaluation =
        crate::commands::select_text(reference, query, context.policy, context.execution).map_err(
            |findings| {
                context.findings(&findings, "reference");
                Exit::of(&findings)
            },
        )?;
    context.findings(&evaluation.warnings, "reference");
    let atoms: Vec<usize> = (&evaluation.selection)
        .into_iter()
        .map(|atom| atom as usize)
        .collect();
    if atoms.is_empty() {
        eprintln!("selection {query:?} matches no atom of the reference");
        return Err(Exit::Consistency);
    }
    Ok(atoms)
}

fn identity(
    atoms: &[usize],
    reference: &molframe::Structure,
    model: &molframe::Structure,
) -> Result<PointMapping, Exit> {
    PointMapping::new(
        atoms.iter().map(|&atom| PointMatch {
            reference: atom,
            model: atom,
        }),
        reference.coordinates().len(),
        model.coordinates().len(),
    )
    .map_err(|error| {
        eprintln!("the model does not share the reference's atoms: {error}");
        Exit::Consistency
    })
}

pub(super) fn region_rmsd(
    metric: MetricChoice,
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Rows {
    let extra = options.extra;
    let (Some(query), Some(alignment)) = (extra.region.as_deref(), extra.alignment) else {
        eprintln!(
            "{} requires --region and --alignment",
            super::command::metric_name(metric)
        );
        return Err(Exit::Usage);
    };
    let atoms = selected(reference, query, context)?;
    let mapping = identity(&atoms, reference, model)?;
    let alignment = match alignment {
        AlignmentChoice::None => ComparisonAlignment::NotRequired,
        AlignmentChoice::Fit => {
            align_mapping(reference.coordinates(), model.coordinates(), &mapping).map_err(
                |error| {
                    eprintln!("region fit failed: {error}");
                    Exit::Consistency
                },
            )?
        }
    };
    let rmsd = match metric {
        MetricChoice::InterfaceRmsd => {
            molframe::compare::interface_rmsd(
                reference.coordinates(),
                model.coordinates(),
                &mapping,
                alignment,
            )
            .0
        }
        _ => {
            molframe::compare::pocket_rmsd(
                reference.coordinates(),
                model.coordinates(),
                &mapping,
                alignment,
            )
            .0
        }
    };
    Ok(vec![(super::command::metric_name(metric), rmsd)])
}

pub(super) fn ligand_rmsd(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Rows {
    use molframe::chemistry::ComponentProvider as _;
    let extra = options.extra;
    let (Some(query), Some(component_id), Some(limit)) = (
        extra.ligand_selection.as_deref(),
        extra.ligand_component.as_deref(),
        extra.ligand_automorphism_limit,
    ) else {
        eprintln!(
            "ligand-rmsd requires --ligand-selection, --ligand-component and \
             --ligand-automorphism-limit"
        );
        return Err(Exit::Usage);
    };
    let source = crate::chemistry::resolve_ccd(options.chemistry.0, options.chemistry.1, context)?;
    let provider = crate::chemistry::load_ccd(source.path, source.version, source.context)?;
    let component = match provider.get(component_id) {
        Ok(Some(component)) => component,
        Ok(None) => {
            eprintln!("component {component_id:?} is absent from the dictionary");
            return Err(Exit::Indeterminate);
        }
        Err(finding) => {
            context.findings(&[finding], "ccd");
            return Err(Exit::Consistency);
        }
    };
    let atoms = selected(reference, query, context)?;
    let moved = match extra.align_on.as_deref() {
        None => model.coordinates().to_vec(),
        Some(fit) => {
            let anchor = selected(reference, fit, context)?;
            let mapping = identity(&anchor, reference, model)?;
            let ComparisonAlignment::Rigid(transform) =
                align_mapping(reference.coordinates(), model.coordinates(), &mapping).map_err(
                    |error| {
                        eprintln!("fit failed: {error}");
                        Exit::Consistency
                    },
                )?
            else {
                return Err(Exit::Consistency);
            };
            model
                .coordinates()
                .iter()
                .map(|position| transform.apply(*position))
                .collect()
        }
    };
    let named = |positions: &[[f32; 3]]| -> Vec<(String, [f32; 3])> {
        atoms
            .iter()
            .filter_map(|&index| {
                let atom = reference
                    .engine()
                    .data()
                    .atom(molframe_core::AtomIndex::new(u32::try_from(index).ok()?))?;
                Some((atom.name()?.to_owned(), *positions.get(index)?))
            })
            .collect()
    };
    let (fixed, mobile) = (named(reference.coordinates()), named(&moved));
    let result = molframe::compare::named_ligand_rmsd(
        &borrowed(&fixed),
        &borrowed(&mobile),
        &component,
        limit,
    )
    .map_err(|error| {
        eprintln!("ligand rmsd failed: {error}");
        Exit::Consistency
    })?;
    Ok(vec![("ligand_rmsd", result.rmsd)])
}

fn borrowed(atoms: &[(String, [f32; 3])]) -> Vec<(&str, [f32; 3])> {
    atoms
        .iter()
        .map(|(name, position)| (name.as_str(), *position))
        .collect()
}
