use super::{DsspError, DsspOptions, secondary_structure};
use molframe_bench::{Sample, structure};
use molframe_chem::{PolymerAtomRole, assign_secondary_structure};
use molframe_core::Presence;
use molframe_core::SecondaryStructure as Ss;
use molframe_core::annotation::{AnnotationColumn, AtomAnnotation};
use molframe_core::structure::Structure;

fn with_roles(structure: &Structure) -> Structure {
    let entries: Vec<_> = structure
        .data()
        .atoms()
        .map(|atom| {
            let role = match atom.name() {
                Some("N") => PolymerAtomRole::PROTEIN_NITROGEN,
                Some("CA") => PolymerAtomRole::PROTEIN_ALPHA_CARBON,
                Some("C") => PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
                Some("O") => PolymerAtomRole::PROTEIN_CARBONYL_OXYGEN,
                _ => PolymerAtomRole::UNKNOWN,
            };
            (role.code(), Presence::Present)
        })
        .collect();
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(AnnotationColumn::from_entries(entries).expect("small column")),
    );
    Structure::new(data)
}

#[test]
fn automatic_and_role_backed_analysis_share_every_exact_state() {
    for sample in [Sample::Tiny, Sample::Small, Sample::Medium, Sample::Large] {
        let source = structure(sample);
        let source = with_roles(&source);
        let automatic = assign_secondary_structure(&source);
        let explicit = secondary_structure(&source, &DsspOptions::default())
            .expect("valid annotated structure");
        assert_eq!(explicit.len(), automatic.len());
        for (residue, kind) in explicit.placements() {
            assert_eq!(
                kind,
                automatic[residue.as_usize()].state,
                "{sample:?} {residue:?}"
            );
        }
    }
}

#[test]
fn proline_and_chain_starts_are_evaluable_without_donating() {
    let source = with_roles(&structure(Sample::Small));
    let table = secondary_structure(&source, &DsspOptions::default()).expect("annotated ubiquitin");
    assert_ne!(table.kind()[0], Ss::Unknown);
    for chain in source.data().chains() {
        for residue in chain.residues().filter(|r| r.name() == Some("PRO")) {
            assert_ne!(table.kind()[residue.index().as_usize()], Ss::Unknown);
        }
    }
}

#[test]
fn explicit_analysis_rejects_missing_annotations_and_invalid_options() {
    let source = structure(Sample::Tiny);
    assert_eq!(
        secondary_structure(&source, &DsspOptions::default()).err(),
        Some(DsspError::MissingRoleAnnotation)
    );
    let invalid = DsspOptions {
        electrostatic_prefactor: f64::NAN,
        ..DsspOptions::default()
    };
    assert_eq!(
        secondary_structure(&source, &invalid).err(),
        Some(DsspError::InvalidOptions)
    );
}

#[test]
fn ambiguous_polymer_roles_are_reported_instead_of_silently_chosen() {
    let source = structure(Sample::Tiny);
    let mut data = source.data().clone();
    let column = AnnotationColumn::from_entries(vec![
        (
            PolymerAtomRole::PROTEIN_NITROGEN.code(),
            Presence::Present
        );
        source.atom_count() as usize
    ])
    .expect("role column");
    data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(column),
    );
    assert!(matches!(
        secondary_structure(&Structure::new(data), &DsspOptions::default()),
        Err(DsspError::AmbiguousRole { .. })
    ));
}
