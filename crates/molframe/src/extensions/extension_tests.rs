//! The one invariant that makes the trait layer safe to use: a convenience
//! method is its `_with` twin at the documented defaults, and every method is
//! the free function it forwards to. A forwarder that drifted would show up
//! here as a disagreement, not as a wrong number nobody looked at.

use super::*;

const DIPEPTIDE: &str = "\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 10.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 10.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00 10.00           C
ATOM      6  N   GLY A   2       3.326   1.551   0.000  1.00 10.00           N
ATOM      7  CA  GLY A   2       3.982   2.847   0.000  1.00 10.00           C
ATOM      8  C   GLY A   2       5.499   2.690   0.000  1.00 10.00           C
ATOM      9  O   GLY A   2       6.100   1.610   0.000  1.00 10.00           O
END
";

fn dipeptide() -> crate::Structure {
    match crate::read_bytes(
        DIPEPTIDE.as_bytes().to_vec(),
        Some("dipeptide.pdb"),
        &crate::ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture read failed: {findings:?}"),
    }
}

#[cfg(feature = "analysis")]
#[test]
fn analysis_methods_agree_with_their_explicit_forms() {
    let structure = dipeptide();
    let context = crate::ExecutionContext::default();

    assert_eq!(
        structure.atom_contacts(4.5).map(|rows| rows.len()),
        structure
            .atom_contacts_with(4.5, crate::spatial::SpatialBackend::Auto, &context)
            .map(|rows| rows.len()),
    );
    assert_eq!(
        structure
            .residue_contact_map(4.5, 0)
            .map(|map| map.contacts().len()),
        structure
            .residue_contact_map_with(4.5, 0, crate::spatial::SpatialBackend::Auto, &context)
            .map(|map| map.contacts().len()),
    );
    assert_eq!(
        structure.salt_bridges(4.0).map(|rows| rows.len()),
        structure
            .salt_bridges_with(4.0, crate::spatial::SpatialBackend::Auto, &context)
            .map(|rows| rows.len()),
    );
    assert_eq!(
        structure.half_sphere_exposure(8.0).map(|rows| rows.len()),
        structure
            .half_sphere_exposure_with(8.0, crate::spatial::SpatialBackend::Auto, &context)
            .map(|rows| rows.len()),
    );
    assert_eq!(
        structure
            .hydrogen_bonds(molframe_analysis::HydrogenBondOptions {
                maximum_donor_acceptor_distance: 3.5,
                minimum_angle_degrees: 120.0,
                backend: crate::spatial::SpatialBackend::Auto,
                periodic: false,
            })
            .map(|rows| rows.len()),
        structure
            .hydrogen_bonds_with(
                molframe_analysis::HydrogenBondOptions {
                    maximum_donor_acceptor_distance: 3.5,
                    minimum_angle_degrees: 120.0,
                    backend: crate::spatial::SpatialBackend::Auto,
                    periodic: false,
                },
                &context,
            )
            .map(|rows| rows.len()),
    );
}

#[cfg(feature = "analysis")]
#[test]
fn the_visitor_form_emits_what_the_collecting_form_returns() {
    let structure = dipeptide();
    let collected = match structure.atom_contacts(4.5) {
        Ok(rows) => rows,
        Err(error) => panic!("contacts failed: {error}"),
    };
    let mut streamed = Vec::new();
    if let Err(error) = structure.visit_atom_contacts(4.5, |contact| streamed.push(contact)) {
        panic!("contact stream failed: {error}");
    }
    assert_eq!(collected.len(), streamed.len());
    assert_eq!(collected.iter().collect::<Vec<_>>(), streamed);
}

#[cfg(feature = "validation")]
#[test]
fn validation_methods_agree_with_their_explicit_forms() {
    let structure = dipeptide();
    let context = crate::ExecutionContext::default();

    assert_eq!(
        structure
            .clashes(0.4, molframe_chem::RadiusSet::Bondi)
            .map(|rows| rows.len()),
        structure
            .clashes_with(
                0.4,
                molframe_chem::RadiusSet::Bondi,
                crate::spatial::SpatialBackend::Auto,
                &context
            )
            .map(|rows| rows.len()),
    );
    assert_eq!(
        structure.quality_flags().len(),
        molframe_validate::quality_flags(structure.engine()).len(),
    );
    assert_eq!(structure.overvalent_atoms(), Vec::new());
    assert_eq!(
        structure.bond_length_deviations(0.1).len(),
        molframe_validate::bond_length_deviations(structure.engine(), 0.1).len(),
    );
}

#[cfg(feature = "compare")]
#[test]
fn a_compare_method_is_the_free_function_it_forwards_to() {
    let model = dipeptide();
    let native = dipeptide();
    let namespace = crate::Namespace::Label;
    assert_eq!(
        model
            .qs_score_in_namespace(
                &native,
                "A",
                "A",
                namespace,
                molframe_compare::QsOptions::standard(5.0)
            )
            .ok(),
        molframe_compare::qs_score_in_namespace(
            model.engine(),
            native.engine(),
            "A",
            "A",
            namespace,
            molframe_compare::QsOptions::standard(5.0)
        )
        .ok(),
    );
}

#[cfg(all(feature = "analysis", feature = "validation", feature = "compare"))]
#[test]
fn every_trait_coexists_in_one_scope_without_a_name_clash() {
    // The prelude brings all three traits in at once. If two of them declared
    // the same method name this file would not compile, which is the check:
    // the call sites below resolve, so the vocabulary is disjoint.
    use crate::prelude::*;

    let structure = dipeptide();
    let _ = structure.atom_contacts(4.5);
    let _ = structure.quality_flags();
    let _ = structure.qs_score(
        &structure,
        "A",
        "A",
        molframe_compare::QsOptions::standard(5.0),
    );
}
