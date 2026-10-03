use crate::bond::{BondRecord, BondTableBuilder};
use crate::diagnostic::Code;
use crate::index::{AtomIndex, ModelIndex, ResidueIndex};
use crate::selection::AtomSelection;
use crate::structure::{CoordinateStore, Structure, StructureData, fixture, validate};
use crate::topology::ModelTable;
use crate::{BondOrder, BondProvenance, SecondarySource, SecondaryStructure};
use std::sync::Arc;

fn annotated() -> Structure {
    let mut data = fixture::sample().data().clone();
    data.secondary_structure = vec![
        SecondaryStructure::AlphaHelix,
        SecondaryStructure::ThreeTenHelix,
        SecondaryStructure::PiHelix,
        SecondaryStructure::OtherHelix,
        SecondaryStructure::PolyProline,
        SecondaryStructure::Unknown,
    ]
    .into();
    data.secondary_source = vec![
        SecondarySource::File,
        SecondarySource::Dssp,
        SecondarySource::CaOnly,
        SecondarySource::File,
        SecondarySource::Dssp,
        SecondarySource::None,
    ]
    .into();
    Structure::new(data)
}

fn select(source: &Structure, atoms: Vec<u32>) -> Structure {
    source
        .materialize(&AtomSelection::from_sorted(atoms))
        .expect("valid materialization")
}

#[test]
fn sparse_materialization_keeps_secondary_states_and_sources_on_their_residues() {
    let source = annotated();
    let selected = select(&source, vec![1, 8, 17]);
    assert_eq!(selected.residue_count(), 3);
    assert_eq!(
        selected.secondary_structure(),
        &[
            SecondaryStructure::AlphaHelix,
            SecondaryStructure::PiHelix,
            SecondaryStructure::PolyProline,
        ]
    );
    assert_eq!(
        selected.secondary_source(),
        &[
            SecondarySource::File,
            SecondarySource::CaOnly,
            SecondarySource::Dssp,
        ]
    );
    for (new, old) in [0, 2, 4].into_iter().enumerate() {
        assert_eq!(
            selected
                .data()
                .topology
                .residues
                .auth_seq_id(ResidueIndex::new(
                    u32::try_from(new).expect("fixture index")
                )),
            source
                .data()
                .topology
                .residues
                .auth_seq_id(ResidueIndex::new(old))
        );
    }
    assert_eq!(source.secondary_structure().len(), 6);
    assert_eq!(
        source.secondary_structure()[4],
        SecondaryStructure::PolyProline
    );
    assert!(validate(selected.data()).is_empty());
}

#[test]
fn selecting_each_residue_preserves_its_secondary_state_and_source() {
    let source = annotated();
    for residue in 0..6 {
        let selected = select(&source, vec![residue * 4 + 2]);
        let row = residue as usize;
        assert_eq!(
            selected.secondary_structure(),
            &source.secondary_structure()[row..=row]
        );
        assert_eq!(
            selected.secondary_source(),
            &source.secondary_source()[row..=row]
        );
    }
}

#[test]
fn all_secondary_classes_survive_remapping_without_reclassification() {
    let source = annotated();
    for state in [
        SecondaryStructure::Unknown,
        SecondaryStructure::Coil,
        SecondaryStructure::AlphaHelix,
        SecondaryStructure::ThreeTenHelix,
        SecondaryStructure::PiHelix,
        SecondaryStructure::OtherHelix,
        SecondaryStructure::BetaBridge,
        SecondaryStructure::Strand,
        SecondaryStructure::Turn,
        SecondaryStructure::Bend,
        SecondaryStructure::PolyProline,
    ] {
        let mut data = source.data().clone();
        data.secondary_structure = vec![state; source.residue_count()].into();
        let selected = select(&Structure::new(data), vec![17]);
        assert_eq!(selected.secondary_structure(), &[state]);
        assert_eq!(selected.secondary_source(), &[SecondarySource::Dssp]);
    }
}

#[test]
fn empty_materialization_discards_all_secondary_rows_without_mutating_the_source() {
    let source = annotated();
    let selected = select(&source, vec![]);
    assert_eq!(selected.atom_count(), 0);
    assert_eq!(selected.residue_count(), 0);
    assert!(selected.secondary_structure().is_empty());
    assert!(selected.secondary_source().is_empty());
    assert_eq!(source.secondary_source().len(), 6);
    let empty = select(&Structure::new(StructureData::empty()), vec![]);
    assert!(empty.secondary_structure().is_empty());
    assert!(validate(empty.data()).is_empty());
}

#[test]
fn unavailable_secondary_columns_stay_unavailable_for_partial_and_empty_selections() {
    let source = fixture::sample();
    for atoms in [vec![4, 20], vec![]] {
        let selected = select(&source, atoms);
        assert!(selected.secondary_structure().is_empty());
        assert!(selected.secondary_source().is_empty());
    }
}

#[test]
fn unchanged_residue_axes_share_secondary_storage_even_when_atoms_are_removed() {
    let source = annotated();
    for selection in [
        AtomSelection::All(source.atom_count()),
        AtomSelection::from_sorted(vec![0, 4, 8, 12, 16, 20]),
    ] {
        let selected = source
            .materialize(&selection)
            .expect("valid materialization");
        assert!(Arc::ptr_eq(
            &selected.data().secondary_structure,
            &source.data().secondary_structure
        ));
        assert!(Arc::ptr_eq(
            &selected.data().secondary_source,
            &source.data().secondary_source
        ));
    }
}

#[test]
fn malformed_secondary_columns_are_rejected_even_if_the_selection_is_empty() {
    let source = annotated();
    for (states, sources) in [(5, 6), (6, 5), (5, 5), (7, 7), (0, 6)] {
        let mut data = source.data().clone();
        data.secondary_structure = vec![SecondaryStructure::Coil; states].into();
        data.secondary_source = vec![SecondarySource::File; sources].into();
        let invalid = Structure::new(data);
        for atoms in [vec![20], vec![]] {
            let result = invalid.materialize(&AtomSelection::from_sorted(atoms));
            assert!(
                matches!(result, Err(findings) if findings.iter().any(|finding| finding.code() == Code::E3011))
            );
        }
    }
}

#[test]
fn dense_model_materialization_compacts_shared_secondary_rows_only_once() {
    let source = annotated();
    let mut data = source.data().clone();
    data.topology.models = ModelTable::default();
    data.topology.models.push(10, 0..2).expect("model 10");
    data.topology.models.push(20, 0..2).expect("model 20");
    let first = data
        .coords
        .block(ModelIndex::new(0))
        .expect("coordinates")
        .clone();
    let second = first
        .as_slice()
        .iter()
        .map(|position| [position[0] + 10.0, position[1], position[2]])
        .collect();
    data.coords = CoordinateStore::Dense {
        frames: vec![first, second],
    };
    let source = Structure::new(data);
    assert!(validate(source.data()).is_empty());
    for atoms in [vec![1, 17], vec![]] {
        let selected = select(&source, atoms.clone());
        assert_eq!(selected.model_count(), 2);
        assert_eq!(selected.residue_count(), atoms.len());
        assert_eq!(selected.secondary_structure().len(), atoms.len());
        assert_eq!(
            selected
                .data()
                .topology
                .models
                .model_num(ModelIndex::new(1)),
            Some(20)
        );
        for model in [0, 1] {
            assert_eq!(
                selected
                    .data()
                    .topology
                    .models
                    .chains(ModelIndex::new(model)),
                Some(0..u32::try_from(selected.chain_count()).expect("fixture chain count"))
            );
            let positions = selected
                .data()
                .coords
                .block(ModelIndex::new(model))
                .expect("retained frame")
                .as_slice();
            for (new, old) in atoms.iter().enumerate() {
                assert_eq!(
                    positions[new].map(f32::to_bits),
                    source
                        .data()
                        .coords
                        .block(ModelIndex::new(model))
                        .expect("source frame")
                        .as_slice()[*old as usize]
                        .map(f32::to_bits)
                );
            }
        }
        if !atoms.is_empty() {
            assert_eq!(
                selected.secondary_structure(),
                &[
                    SecondaryStructure::AlphaHelix,
                    SecondaryStructure::PolyProline
                ]
            );
            assert_eq!(
                selected.secondary_source(),
                &[SecondarySource::File, SecondarySource::Dssp]
            );
        }
    }
}

#[test]
fn model_local_residues_keep_secondary_identity_when_an_entire_model_loses_its_atoms() {
    let source = annotated();
    let mut data = source.data().clone();
    data.topology.models = ModelTable::default();
    data.topology.models.push(10, 0..1).expect("model 10");
    data.topology.models.push(20, 1..2).expect("model 20");
    let source = Structure::new(data);
    assert!(validate(source.data()).is_empty());
    for atoms in [vec![1, 17], vec![17]] {
        let selected = select(&source, atoms.clone());
        let expected: Vec<_> = atoms
            .iter()
            .map(|atom| source.secondary_structure()[*atom as usize / 4])
            .collect();
        assert_eq!(selected.secondary_structure(), expected);
        assert_eq!(
            selected
                .data()
                .topology
                .models
                .model_num(ModelIndex::new(1)),
            Some(20)
        );
        let expected_first = if atoms.len() == 2 { 0..1 } else { 0..0 };
        assert_eq!(
            selected.data().topology.models.chains(ModelIndex::new(0)),
            Some(expected_first)
        );
    }
}

#[test]
fn materialization_remaps_retained_bond_endpoints_order_and_provenance() {
    let source = annotated();
    let mut data = source.data().clone();
    let mut bonds = BondTableBuilder::new();
    for (atom_a, atom_b, order, provenance) in [
        (0, 17, BondOrder::Single, BondProvenance::File),
        (4, 17, BondOrder::Double, BondProvenance::File),
        (17, 20, BondOrder::Single, BondProvenance::User),
    ] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(atom_a),
            atom_b: AtomIndex::new(atom_b),
            order,
            provenance,
        });
    }
    data.bonds = bonds.finish();
    let source = Structure::new(data);
    let selected = select(&source, vec![4, 17, 20]);
    let retained: Vec<_> = selected
        .data()
        .bonds
        .iter()
        .map(|bond| {
            (
                bond.atom_a.get(),
                bond.atom_b.get(),
                bond.order,
                bond.provenance,
            )
        })
        .collect();
    assert_eq!(
        retained,
        vec![
            (0, 1, BondOrder::Double, BondProvenance::File),
            (1, 2, BondOrder::Single, BondProvenance::User)
        ]
    );
    assert_eq!(source.data().bonds.len(), 3);
    assert!(select(&source, vec![]).data().bonds.is_available());
    assert!(
        !select(&fixture::sample(), vec![4])
            .data()
            .bonds
            .is_available()
    );
}

#[test]
fn ragged_materialization_is_rejected_without_changing_child_secondary_rows() {
    let child = annotated();
    let mut data = StructureData::empty();
    data.coords = CoordinateStore::Ragged {
        models: vec![child.clone()],
    };
    let source = Structure::new(data);
    assert!(
        matches!(source.materialize(&AtomSelection::All(0)), Err(findings) if findings.iter().any(|finding| finding.code() == Code::E6003))
    );
    assert_eq!(
        child.secondary_structure()[4],
        SecondaryStructure::PolyProline
    );
}
