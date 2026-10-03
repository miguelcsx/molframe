use super::*;
use crate::{ComponentAtom, ComponentBond, ComponentKind, MemoryProvider};
use molframe_core::annotation::AnnotationColumn;
use molframe_core::bond::{BondOrder, BondTableBuilder};
use molframe_core::io::{InputBuffer, ReadOptions};

fn component() -> Component {
    Component {
        id: "LIG".into(),
        name: "methanol".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: [
            ("C", molframe_core::Element::CARBON),
            ("O", molframe_core::Element::OXYGEN),
            ("H", molframe_core::Element::HYDROGEN),
        ]
        .into_iter()
        .map(|(name, element)| ComponentAtom {
            name: name.into(),
            alternate_name: None,
            element,
            charge: 0,
            aromatic: false,
            leaving: false,
            stereo: None,
        })
        .collect::<Vec<_>>()
        .into(),
        bonds: [("C", "O"), ("O", "H")]
            .into_iter()
            .map(|(a, b)| ComponentBond {
                atom_a: a.into(),
                atom_b: b.into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            })
            .collect::<Vec<_>>()
            .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn provider() -> MemoryProvider {
    MemoryProvider::new(DictionaryVersion::new("fixture"), [component()]).unwrap()
}

fn structure() -> Structure {
    let text = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 O O LIG A 1 1.4 0 0\nHETATM 2 C C LIG A 1 0 0 0\n";
    let parsed = molframe_cif::read(
        &InputBuffer::from_bytes(text.as_bytes().to_vec()),
        &ReadOptions::new(),
    )
    .unwrap()
    .0;
    let mut data = parsed.data().clone();
    data.bonds = BondTableBuilder::new().finish();
    Structure::new(data)
}

#[test]
fn missing_hydrogen_charge_is_projected_in_observed_atom_order() {
    let expected = crate::component_peoe_charges(&component(), PeoeOptions::default()).unwrap();
    let result = partial_charges(&structure(), &provider(), PeoeOptions::default()).unwrap();
    assert!((result.values[0] - expected[1] - expected[2]).abs() < 1.0e-12);
    assert!((result.values[1] - expected[0]).abs() < 1.0e-12);
    assert!(result.values.iter().sum::<f64>().abs() < 1.0e-12);
    assert_eq!(
        result.source,
        ChargeSource::Peoe {
            dictionary: DictionaryVersion::new("fixture"),
            options: PeoeOptions::default()
        }
    );
}

#[test]
fn a_complete_file_column_wins_even_with_an_unknown_dictionary_component() {
    let original = structure();
    let mut data = original.data().clone();
    data.annotations.insert(
        PARTIAL_CHARGE_ANNOTATION,
        AtomAnnotation::Real(
            AnnotationColumn::from_entries([
                (0.125, Presence::Present),
                (-0.25, Presence::Present),
            ])
            .unwrap(),
        ),
    );
    let unknown = MemoryProvider::new(DictionaryVersion::new("empty"), []).unwrap();
    let result = partial_charges(&Structure::new(data), &unknown, PeoeOptions::default()).unwrap();
    assert_eq!(
        result,
        PartialCharges {
            values: vec![0.125, -0.25],
            source: ChargeSource::File
        }
    );
}

#[test]
fn partial_file_columns_are_not_mixed_with_calculated_charges() {
    let original = structure();
    let mut data = original.data().clone();
    data.annotations.insert(
        PARTIAL_CHARGE_ANNOTATION,
        AtomAnnotation::Real(
            AnnotationColumn::from_entries([(0.125, Presence::Present), (0.0, Presence::Unknown)])
                .unwrap(),
        ),
    );
    assert_eq!(
        partial_charges(&Structure::new(data), &provider(), PeoeOptions::default()),
        Err(PartialChargeError::InvalidFileCharges)
    );
}

#[test]
fn unresolved_components_and_absent_heavy_atoms_fail_explicitly() {
    let empty = MemoryProvider::new(DictionaryVersion::new("empty"), []).unwrap();
    assert_eq!(
        partial_charges(&structure(), &empty, PeoeOptions::default()),
        Err(PartialChargeError::UnknownComponent("LIG".into()))
    );
    let mut component = component();
    let mut atoms = component.atoms.to_vec();
    let mut absent = atoms[0].clone();
    absent.name = "C2".into();
    atoms.push(absent);
    component.atoms = atoms.into();
    let incomplete = MemoryProvider::new(DictionaryVersion::new("fixture"), [component]).unwrap();
    assert_eq!(
        partial_charges(&structure(), &incomplete, PeoeOptions::default()),
        Err(PartialChargeError::MissingAtom {
            component: "LIG".into(),
            atom: "C2".into()
        })
    );
}

#[test]
fn complete_4hhb_protein_charges_resolve_linked_ccd_chemistry_in_stable_order() {
    let source = molframe_bench::structure_from_pdb(molframe_bench::Sample::Medium).unwrap();
    let (provider, _) = crate::read_ccd(
        &InputBuffer::from_bytes(include_bytes!("../../data/CCD-amino-acids.cif").to_vec()),
        DictionaryVersion::new("wwPDB-2026-10-03"),
    )
    .unwrap();
    let selected = molframe_core::selection::AtomSelection::from_sorted(
        source
            .data()
            .atoms()
            .filter(|atom| {
                atom.residue()
                    .and_then(ResidueRef::name)
                    .is_some_and(|name| provider.get(name).unwrap().is_some())
            })
            .map(|atom| atom.index().get())
            .collect(),
    );
    let protein = source.materialize(&selected).unwrap();
    let protein = crate::perceive_bonds(&protein).unwrap();
    let result = partial_charges(&protein, &provider, PeoeOptions::default()).unwrap();
    assert_eq!(result.values.len(), protein.atom_count() as usize);
    assert!(result.values.iter().all(|charge| charge.is_finite()));
    let lysine_nz: Vec<_> = protein
        .data()
        .atoms()
        .filter(|atom| {
            atom.name() == Some("NZ") && atom.residue().and_then(ResidueRef::name) == Some("LYS")
        })
        .collect();
    assert!(!lysine_nz.is_empty());
    for atom in lysine_nz {
        assert!(result.values[atom.index().as_usize()] > 0.0);
    }
}
