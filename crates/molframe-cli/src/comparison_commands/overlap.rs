//! Scores of how much of a reference's contacts a model reproduces.

use super::command::ComparisonOptions;
use crate::exit::Exit;
use crate::report::Context;
use molframe::spatial::SpatialBackend;

type Rows = Result<Vec<(&'static str, f64)>, Exit>;

fn missing(metric: &str, flags: &str) -> Exit {
    eprintln!("{metric} requires {flags}");
    Exit::Usage
}

fn failed(metric: &str, error: &dyn std::fmt::Display) -> Exit {
    eprintln!("{metric} failed: {error}");
    Exit::Consistency
}

pub(super) fn qs(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Rows {
    let (Some(first), Some(second), Some(cutoff)) =
        (options.receptor, options.ligand, options.contact_distance)
    else {
        return Err(missing("qs", "--receptor, --ligand and --contact-distance"));
    };
    let qs_options = molframe::compare::QsOptions::standard(cutoff);
    let score = match options.mapped {
        None => molframe::compare::qs_score_in_namespace(
            model.engine(),
            reference.engine(),
            first,
            second,
            context.policy.identifiers,
            qs_options,
        )
        .map_err(|error| failed("qs", &error))?,
        Some(mapped) => {
            use molframe::CompareExt as _;
            let provider = crate::chemistry::load_ccd(mapped.ccd, mapped.ccd_version, context)?;
            let mapping = molframe::compare::MappingOptions {
                provider: &provider,
                namespace: context.policy.identifiers,
                scoring: mapped.scoring,
                min_identity: mapped.min_identity,
                automorphism_limit: mapped.automorphism_limit,
            };
            model
                .mapped_qs_score(reference, first, second, &mapping, qs_options)
                .map_err(|error| failed("qs", &error))?
                .0
        }
    };
    Ok(vec![("qs", score)])
}

fn residue_pairs(
    structure: &molframe::Structure,
    cutoff: f32,
    min_separation: u32,
    context: Context,
) -> Result<Vec<(u32, u32)>, Exit> {
    let map = molframe::analysis::residue_contact_map(
        structure.engine(),
        cutoff,
        min_separation,
        SpatialBackend::Auto,
        context.execution,
    )
    .map_err(|error| failed("contact similarity", &error))?;
    Ok(map
        .contacts()
        .iter()
        .map(|contact| (contact.first.get(), contact.second.get()))
        .collect())
}

pub(super) fn contact_similarity(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Rows {
    let (Some(cutoff), Some(separation)) = (options.contact_distance, options.extra.min_separation)
    else {
        return Err(missing(
            "contact-similarity",
            "--contact-distance and --min-separation",
        ));
    };
    let from_model = residue_pairs(model, cutoff, separation, context)?;
    let from_reference = residue_pairs(reference, cutoff, separation, context)?;
    let similarity = molframe::compare::contact_map_similarity(&from_reference, &from_model);
    Ok(vec![
        ("contact_similarity", similarity.jaccard),
        ("contact_similarity_shared", count(similarity.shared)),
        ("contact_similarity_union", count(similarity.union)),
    ])
}

#[allow(clippy::cast_precision_loss)]
const fn count(value: usize) -> f64 {
    value as f64
}

/// Per-atom radii and residue identifiers, which the contact-area surface needs.
fn surface_inputs(
    structure: &molframe::Structure,
    radii: molframe::chemistry::RadiusSet,
) -> Result<(Vec<f32>, Vec<u32>), Exit> {
    let mut radius = Vec::new();
    let mut residue = Vec::new();
    for atom in structure.engine().data().atoms() {
        let (Some(element), Some(owner)) = (atom.element(), atom.residue()) else {
            eprintln!("cad needs an element and a residue for every atom");
            return Err(Exit::Consistency);
        };
        let Some(value) = molframe::chemistry::vdw_radius(element, radii) else {
            eprintln!("cad: the radius set has no radius for {}", element.symbol());
            return Err(Exit::Consistency);
        };
        radius.push(value);
        residue.push(owner.index().get());
    }
    Ok((radius, residue))
}

pub(super) fn cad(
    model: &molframe::Structure,
    reference: &molframe::Structure,
    options: &ComparisonOptions<'_>,
    context: Context,
) -> Rows {
    let extra = options.extra;
    let (Some(probe), Some(density), Some(radii)) =
        (extra.cad_probe, extra.cad_density, extra.radii)
    else {
        return Err(missing("cad", "--cad-probe, --cad-density and --radii"));
    };
    let areas = |structure: &molframe::Structure| {
        let (radius, residue) = surface_inputs(structure, radii.into())?;
        molframe::compare::cad_contact_areas(
            structure.coordinates(),
            &radius,
            &residue,
            probe,
            density,
            context.execution,
        )
        .map_err(|error| failed("cad", &error))
    };
    let score = molframe::compare::cad_score(&areas(reference)?, &areas(model)?)
        .map_err(|error| failed("cad", &error))?;
    Ok(vec![
        ("cad", score.score),
        ("cad_reference_area", score.reference_area),
        ("cad_lost_area", score.lost_area),
    ])
}
