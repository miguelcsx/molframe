use super::{SurfaceSampling, definition_contacts_kernel};
use crate::analyse_structure;
use molframe_cif::read;
use molframe_core::contract::{
    AnalysisPolicy, ContactDefinition, Indeterminacy, MissingPolicy, PolicyField, RadiiSet,
};
use molframe_core::{ExecutionContext, InputBuffer, ReadOptions, Structure};
use molframe_spatial::SpatialBackend;
use std::fmt::Write;

const HEADER: &str = "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn carbons(rows: &[(&str, f32)]) -> Structure {
    let mut text = HEADER.to_owned();
    for (index, (symbol, x)) in rows.iter().enumerate() {
        let number = index + 1;
        let _ = writeln!(
            text,
            "ATOM {number} {symbol} C{number} GLY A {number} {x} 0 0"
        );
    }
    let input = InputBuffer::from_bytes(text.into_bytes());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn policy(definition: ContactDefinition, radii: RadiiSet) -> AnalysisPolicy {
    AnalysisPolicy {
        contact_def: definition,
        vdw_radii: radii,
        ..AnalysisPolicy::default()
    }
}

fn count(structure: &Structure, policy: &AnalysisPolicy) -> usize {
    let kernel = definition_contacts_kernel(
        SpatialBackend::BruteForce,
        SurfaceSampling {
            tolerance: 0.2,
            density: 4.0,
            minimum_area: 0.25,
        },
    );
    match analyse_structure(structure, policy, &kernel, &ExecutionContext::default()) {
        Ok(result) => result.value().map_or(0, crate::ContactTable::len),
        Err(error) => panic!("analysis failed: {error}"),
    }
}

fn carbon_radius(set: RadiiSet) -> f32 {
    let Ok(chem_set) = set.name().parse::<molframe_chem::RadiusSet>() else {
        panic!("the radius set is known to chemistry");
    };
    let Some(radius) = molframe_chem::vdw_radius(molframe_core::Element::CARBON, chem_set) else {
        panic!("carbon has a radius in every set");
    };
    radius
}

#[test]
fn the_distance_is_measured_against_the_radii_the_policy_names() {
    for set in [
        RadiiSet::Bondi,
        RadiiSet::AmberUnited,
        RadiiSet::Charmm,
        RadiiSet::Alvarez,
    ] {
        let touching = 2.0 * carbon_radius(set);
        let near = two_carbons(touching - 0.01);
        let far = two_carbons(touching + 0.01);
        let definition = ContactDefinition::DistanceCutoff { tolerance: 0.0 };
        assert_eq!(
            count(&near, &policy(definition, set)),
            1,
            "{set:?}: just inside"
        );
        assert_eq!(
            count(&far, &policy(definition, set)),
            0,
            "{set:?}: just outside"
        );
    }
}

fn two_carbons(separation: f32) -> Structure {
    carbons(&[("C", 0.0), ("C", separation)])
}

#[test]
fn the_tolerance_widens_the_contact_by_exactly_that_much() {
    let structure = two_carbons(2.0 * carbon_radius(RadiiSet::Bondi) + 0.3);
    let at = |tolerance| {
        count(
            &structure,
            &policy(
                ContactDefinition::DistanceCutoff { tolerance },
                RadiiSet::Bondi,
            ),
        )
    };
    assert_eq!(at(0.25), 0);
    assert_eq!(at(0.35), 1);
}

#[test]
fn an_atom_with_no_radius_is_missing_coverage_and_the_policy_decides_what_that_means() {
    let structure = carbons(&[("C", 0.0), ("C", 3.0), ("?", 6.0)]);
    let kernel = definition_contacts_kernel(
        SpatialBackend::BruteForce,
        SurfaceSampling {
            tolerance: 0.2,
            density: 4.0,
            minimum_area: 0.25,
        },
    );
    let context = ExecutionContext::default();
    let reported = AnalysisPolicy {
        missing_atoms: MissingPolicy::Report,
        ..AnalysisPolicy::default()
    };
    let Ok(result) = analyse_structure(&structure, &reported, &kernel, &context) else {
        panic!("a reported shortfall is a result");
    };
    assert_eq!(result.coverage.missing, 1);
    assert_eq!(result.value().map(crate::ContactTable::len), Some(1));

    let refused = AnalysisPolicy {
        missing_atoms: MissingPolicy::Indeterminate,
        ..AnalysisPolicy::default()
    };
    let Ok(result) = analyse_structure(&structure, &refused, &kernel, &context) else {
        panic!("an indeterminate policy is a result");
    };
    assert!(result.value().is_none());
    assert!(matches!(
        result.indeterminacy(),
        Some(Indeterminacy::Frame { .. })
    ));
}

#[test]
fn a_surface_contact_needs_exposed_surface_on_both_atoms() {
    let structure = two_carbons(3.5);
    let surface = policy(
        ContactDefinition::SurfaceBased { probe: 1.4 },
        RadiiSet::Bondi,
    );
    assert_eq!(count(&structure, &surface), 1);
    let far = two_carbons(30.0);
    assert_eq!(count(&far, &surface), 0);
}

#[test]
fn the_kernel_records_the_two_decisions_it_applied() {
    let structure = two_carbons(3.0);
    let kernel = definition_contacts_kernel(
        SpatialBackend::BruteForce,
        SurfaceSampling {
            tolerance: 0.2,
            density: 4.0,
            minimum_area: 0.25,
        },
    );
    let Ok(result) = analyse_structure(
        &structure,
        &AnalysisPolicy::default(),
        &kernel,
        &ExecutionContext::default(),
    ) else {
        panic!("a plain analysis");
    };
    let Some(reads) = result.provenance.policy_reads() else {
        panic!("reads are recorded");
    };
    assert!(reads.contains(&PolicyField::ContactDef) && reads.contains(&PolicyField::VdwRadii));
}
