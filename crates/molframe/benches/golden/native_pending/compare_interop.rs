use std::fmt::Debug;
use std::sync::Arc;

use criterion::{BenchmarkGroup, Throughput, black_box};
use molframe::chemistry::{Component, ComponentKind};
use molframe::{DictionaryVersion, Element, ReadOptions};

trait BenchRequired<T> {
    fn required(self, context: &str) -> T;
}

impl<T, E: Debug> BenchRequired<T> for Result<T, E> {
    fn required(self, context: &str) -> T {
        match self {
            Ok(value) => value,
            Err(error) => panic!("{context}: {error:?}"),
        }
    }
}

impl<T> BenchRequired<T> for Option<T> {
    fn required(self, context: &str) -> T {
        match self {
            Some(value) => value,
            None => panic!("{context}"),
        }
    }
}

pub(super) fn register(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    bench_gw_027(group);
    bench_gw_029(group);
    bench_gw_034(group);
    bench_gw_036(group);
    bench_gw_041(group);
}

fn bench_gw_027(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
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
    group.throughput(Throughput::Elements(4));
    group.bench_function("GW-027", |b| {
        b.iter(|| {
            let assignment = molframe::compare::assign_chains(
                reference.engine(),
                target.engine(),
                &provider,
                molframe::Namespace::Label,
                molframe::sequence::Scoring::simple(),
                1.0,
            )
            .required("GW-027 failed");
            black_box((assignment.primary.len(), assignment.alternatives.len()));
        });
    });
}

fn bench_gw_029(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let component = symmetric_component();
    group.throughput(Throughput::Elements(3));
    group.bench_function("GW-029", |b| {
        b.iter(|| {
            let result = molframe::compare::ligand_symmetry_rmsd(
                &[[0.0, 0.0, 0.0], [-1.2, 0.0, 0.0], [1.2, 0.0, 0.0]],
                &[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [-1.2, 0.0, 0.0]],
                &component,
                4,
            )
            .required("GW-029 failed");
            black_box(result.rmsd);
        });
    });
}

fn bench_gw_034(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
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
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-034", |b| {
        b.iter(|| {
            let tensor = molframe::interop::DlpackTensor::coordinates(structure.engine())
                .required("GW-034 failed");
            black_box((tensor.cost(), tensor.as_managed().is_some()));
        });
    });
}

fn bench_gw_036(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let dataset = molframe::interop::Dataset::new(vec![
        entry("a", "AAAA", "2020-01-01"),
        entry("b", "AAAA", "2020-01-02"),
        entry("c", "GGGG", "2021-01-01"),
        entry("d", "GGGG", "2021-01-02"),
    ])
    .required("GW-036 fixture failed");
    let ratios =
        molframe::interop::SplitRatios::new(0.5, 0.25, 0.25).required("GW-036 ratios failed");
    group.throughput(Throughput::Elements(4));
    group.bench_function("GW-036", |b| {
        b.iter(|| {
            let split = dataset
                .split(&molframe::interop::SplitOptions {
                    strategy: molframe::interop::SplitStrategy::SequenceIdentity { threshold: 1.0 },
                    ratios,
                })
                .required("GW-036 failed");
            black_box((split.train.len(), split.validation.len(), split.test.len()));
        });
    });
}

fn bench_gw_041(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let input = b"golden-provenance-input";
    let policy = molframe::AnalysisPolicy::default();
    let provenance = molframe::Provenance::new(&policy)
        .with_source(molframe::SourceRef::Memory)
        .with_input_fingerprint(molframe_core::contract::Fingerprint::of(input));
    group.throughput(Throughput::Bytes(input.len() as u64));
    group.bench_function("GW-041", |b| {
        b.iter(|| {
            let replay = molframe_core::contract::reexecute_from_provenance(
                &provenance,
                input,
                molframe_core::contract::ReexecutionEnvironment::current(),
                |bytes, replay_policy| (bytes.len(), replay_policy.fingerprint()),
            )
            .required("GW-041 failed");
            black_box(replay.value);
        });
    });
}

fn read(source: &str) -> molframe::Structure {
    molframe::read_bytes(
        source.as_bytes().to_vec(),
        Some("golden.cif"),
        &ReadOptions::new(),
    )
    .required("golden benchmark fixture failed")
    .0
}

fn provider() -> molframe::chemistry::MemoryProvider {
    molframe::chemistry::MemoryProvider::new(
        DictionaryVersion::new("golden-ccd"),
        [component("GLY", b'G'), component("ALA", b'A')],
    )
    .required("CCD provider failed")
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
            molframe::chemistry::ComponentBond {
                atom_a: "C".into(),
                atom_b: "O1".into(),
                order: molframe::BondOrder::Double,
                aromatic: false,
                stereo: None,
            },
            molframe::chemistry::ComponentBond {
                atom_a: "C".into(),
                atom_b: "O2".into(),
                order: molframe::BondOrder::Double,
                aromatic: false,
                stereo: None,
            },
        ]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn atom(name: &str, element: Element) -> molframe::chemistry::ComponentAtom {
    molframe::chemistry::ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}

fn entry(id: &str, sequence: &str, date: &str) -> molframe::interop::ManifestEntry {
    molframe::interop::ManifestEntry {
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
