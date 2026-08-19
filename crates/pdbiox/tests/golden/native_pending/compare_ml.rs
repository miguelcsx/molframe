use std::sync::Arc;

use pdbiox::{Component, ComponentKind, DictionaryVersion, Element, ReadOptions};

#[test]
fn gw_027_maps_homomeric_chains_and_reports_alternatives() {
    let reference = read(
        r"data_reference
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY A 1 0 0 0
2 C CA GLY A 2 1 0 0
3 C CA GLY B 1 0 0 0
4 C CA GLY B 2 1 0 0
",
    );
    let target = read(
        r"data_target
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA GLY X 1 0 0 0
2 C CA GLY X 2 1 0 0
3 C CA GLY Y 1 0 0 0
4 C CA GLY Y 2 1 0 0
",
    );
    let provider = provider();
    let assignment = pdbiox::compare::assign_chains(
        &reference,
        &target,
        &provider,
        pdbiox::Namespace::Label,
        pdbiox::seq::Scoring::simple(),
        1.0,
    )
    .unwrap_or_else(|finding| panic!("chain mapping failed: {finding}"));
    assert_eq!(assignment.primary.len(), 2);
    assert!(
        assignment
            .primary
            .iter()
            .all(|mapping| (mapping.identity - 1.0).abs() < f64::EPSILON)
    );
    assert_eq!(assignment.alternatives.len(), 2);
    assert!(
        assignment
            .alternatives
            .iter()
            .all(|alternative| (alternative.identity - 1.0).abs() < f64::EPSILON)
    );
}

#[test]
fn gw_029_uses_chemical_symmetry_for_ligand_pose_comparison() {
    let component = symmetric_component();
    let result = pdbiox::compare::ligand_symmetry_rmsd(
        &[[0.0, 0.0, 0.0], [-1.2, 0.0, 0.0], [1.2, 0.0, 0.0]],
        &[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [-1.2, 0.0, 0.0]],
        &component,
        4,
    )
    .unwrap_or_else(|error| panic!("ligand symmetry comparison failed: {error}"));
    assert!(result.rmsd < f64::EPSILON);
    assert_eq!(result.mapping.reference_to_model.as_ref(), &[0, 2, 1]);
}

#[test]
fn gw_034_exports_a_mutation_isolated_dlpack_tensor() {
    let structure = read(
        r"data_tensor
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA ALA A 1 1 2 3
ATOM 2 N N ALA A 1 4 5 6
",
    );
    let tensor = pdbiox::DlpackTensor::coordinates(&structure)
        .unwrap_or_else(|error| panic!("DLPack export failed: {error}"));
    let (data, shape, dtype_bits) = {
        let managed = tensor
            .as_managed()
            .unwrap_or_else(|| panic!("DLPack manager missing"));
        let shape = unsafe { std::slice::from_raw_parts(managed.dl_tensor.shape, 2) };
        (
            managed.dl_tensor.data.cast::<f32>(),
            [shape[0], shape[1]],
            managed.dl_tensor.dtype.bits,
        )
    };
    assert_eq!(shape, [2, 3]);
    assert_eq!(dtype_bits, 32);
    assert_eq!(tensor.cost(), pdbiox::ExportCost::Copy);
    assert_ne!(data.cast_const(), structure.positions().as_ptr().cast());
    unsafe { data.write(9.0) };
    assert_eq!(structure.positions()[0][0].to_bits(), 1.0_f32.to_bits());
}

#[test]
fn gw_036_splits_a_manifest_by_sequence_identity_without_loading_coordinates() {
    let dataset = pdbiox::Dataset::new(vec![
        entry("a", "AAAA", "2020-01-01"),
        entry("b", "AAAA", "2020-01-02"),
        entry("c", "GGGG", "2021-01-01"),
        entry("d", "GGGG", "2021-01-02"),
    ])
    .unwrap_or_else(|error| panic!("dataset fixture failed: {error}"));
    let ratios = pdbiox::SplitRatios::new(0.5, 0.25, 0.25)
        .unwrap_or_else(|error| panic!("split ratios failed: {error}"));
    let split = dataset
        .split(&pdbiox::SplitOptions {
            strategy: pdbiox::SplitStrategy::SequenceIdentity { threshold: 1.0 },
            ratios,
        })
        .unwrap_or_else(|error| panic!("sequence split failed: {error}"));
    assert_eq!(
        split.train.len() + split.validation.len() + split.test.len(),
        4
    );
    assert!(
        split
            .train
            .entries()
            .any(|entry| entry.sequence.as_deref() == Some("AAAA"))
    );
    assert!(
        split
            .validation
            .entries()
            .chain(split.test.entries())
            .all(|entry| entry.path.extension().is_none())
    );
}

#[test]
fn gw_041_reexecutes_only_with_matching_provenance() {
    let input = b"golden-provenance-input";
    let policy = pdbiox::AnalysisPolicy::default();
    let provenance = pdbiox::Provenance::new(&policy)
        .with_source(pdbiox::SourceRef::Memory)
        .with_input_fingerprint(pdbiox::core::contract::Fingerprint::of(input));
    let replay = pdbiox::core::contract::reexecute_from_provenance(
        &provenance,
        input,
        pdbiox::core::contract::ReexecutionEnvironment::current(),
        |bytes, replay_policy| (bytes.len(), replay_policy.fingerprint()),
    )
    .unwrap_or_else(|error| panic!("provenance replay failed: {error}"));
    assert_eq!(replay.value.0, input.len());
    assert_eq!(replay.value.1, policy.fingerprint());
    assert_eq!(replay.provenance.fingerprint(), provenance.fingerprint());
    assert!(
        pdbiox::core::contract::reexecute_from_provenance(
            &provenance,
            b"different-input",
            pdbiox::core::contract::ReexecutionEnvironment::current(),
            |_, _| (),
        )
        .is_err()
    );
}

fn read(source: &str) -> pdbiox::Structure {
    match pdbiox::read_bytes(
        source.as_bytes().to_vec(),
        Some("golden.cif"),
        &ReadOptions::new(),
    ) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("golden native fixture failed: {findings:?}"),
    }
}

fn provider() -> pdbiox::MemoryProvider {
    pdbiox::MemoryProvider::new(
        DictionaryVersion::new("golden-ccd"),
        [component("GLY", b'G'), component("ALA", b'A')],
    )
    .unwrap_or_else(|finding| panic!("CCD provider failed: {finding}"))
}

fn component(id: &str, code: u8) -> Component {
    Component {
        id: id.into(),
        name: id.into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(code),
        formula: None,
        atoms: Arc::from([]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn symmetric_component() -> Component {
    Component {
        id: "CO2".into(),
        name: "carbon dioxide".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: Some("CO2".into()),
        atoms: Arc::from([
            atom("C", Element::CARBON),
            atom("O1", Element::OXYGEN),
            atom("O2", Element::OXYGEN),
        ]),
        bonds: Arc::from([
            pdbiox::ComponentBond {
                atom_a: "C".into(),
                atom_b: "O1".into(),
                order: pdbiox::BondOrder::Double,
                aromatic: false,
                stereo: None,
            },
            pdbiox::ComponentBond {
                atom_a: "C".into(),
                atom_b: "O2".into(),
                order: pdbiox::BondOrder::Double,
                aromatic: false,
                stereo: None,
            },
        ]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn atom(name: &str, element: Element) -> pdbiox::ComponentAtom {
    pdbiox::ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}

fn entry(id: &str, sequence: &str, date: &str) -> pdbiox::ManifestEntry {
    pdbiox::ManifestEntry {
        id: id.into(),
        path: id.into(),
        atom_count: 100,
        resolution: None,
        method: None,
        deposition_date: Some(date.into()),
        sequence: Some(sequence.into()),
        structure_cluster: None,
        tags: Vec::new(),
        statistics: std::collections::BTreeMap::new(),
    }
}
