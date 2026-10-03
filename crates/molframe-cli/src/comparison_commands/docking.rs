//! `DockQ`, scored directly or after mapping chains by sequence.

use super::command::{ComparisonOptions, MappedScoring};
use crate::exit::Exit;
use crate::report::Context;

pub(super) fn measure_dockq(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Result<Vec<(&'static str, f64)>, Exit> {
    let Some(receptor) = options.receptor else {
        eprintln!("--receptor is required when dock-q is selected");
        return Err(Exit::Usage);
    };
    let Some(ligand) = options.ligand else {
        eprintln!("--ligand is required when dock-q is selected");
        return Err(Exit::Usage);
    };
    let (Some(contact_distance), Some(ligand_scale), Some(interface_scale)) = (
        options.contact_distance,
        options.ligand_scale,
        options.interface_scale,
    ) else {
        eprintln!(
            "--contact-distance, --ligand-scale and --interface-scale are required when dock-q is selected"
        );
        return Err(Exit::Usage);
    };
    let dockq_options = molframe::compare::DockQOptions {
        contact_distance,
        ligand_scale,
        interface_scale,
    };
    let scored = match options.mapped {
        None => molframe::compare::dockq(
            model.engine(),
            reference.engine(),
            receptor,
            ligand,
            dockq_options,
        )
        .map_err(|error| error.to_string()),
        Some(mapped) => mapped_dockq(
            model,
            reference,
            (receptor, ligand),
            mapped,
            dockq_options,
            context,
        ),
    };
    scored
        .map(|result| {
            vec![
                ("dockq", result.score),
                ("dockq_fnat", result.fnat),
                ("dockq_ligand_rmsd", result.ligand_rmsd),
                ("dockq_interface_rmsd", result.interface_rmsd),
            ]
        })
        .map_err(|error| {
            eprintln!("DockQ failed: {error}");
            Exit::Consistency
        })
}

/// `DockQ` after matching chains, residues and atoms by sequence and name.
fn mapped_dockq(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    (receptor, ligand): (&str, &str),
    mapped: MappedScoring<'_>,
    options: molframe::compare::DockQOptions,
    context: Context,
) -> Result<molframe::compare::DockQ, String> {
    use molframe::CompareExt;
    let provider = crate::chemistry::load_ccd(mapped.ccd, mapped.ccd_version, context)
        .map_err(|_| "the Chemical Component Dictionary could not be loaded".to_owned())?;
    let mapping = molframe::compare::MappingOptions {
        provider: &provider,
        namespace: context.policy.identifiers,
        scoring: mapped.scoring,
        min_identity: mapped.min_identity,
        automorphism_limit: mapped.automorphism_limit,
    };
    model
        .mapped_dockq(reference, receptor, ligand, &mapping, options)
        .map(|(score, _)| score)
        .map_err(|error| error.to_string())
}
