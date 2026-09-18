//! Criterion measurements for the remaining executable native golden workflows.

use std::collections::BTreeMap;

use criterion::{BenchmarkGroup, Throughput, black_box};
use molframe::{
    AltlocPolicy, AnalysisPolicy, AtomSelection, InputBuffer, ModelChoice, ModelCifExt,
    ReadOptions, Rigid,
};
use molframe_bench::{Sample, coordinates, structure};

const UNKNOWN_CATEGORY_CIF: &str = "data_unknown\n\
_custom.note 'keep this category'\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n";

const MODEL_CIF: &str = "data_model\n\
loop_\n_ma_qa_metric.id\n_ma_qa_metric.name\n_ma_qa_metric.type\n_ma_qa_metric.mode\n\
1 pLDDT pLDDT local\n2 PAE PAE local-pairwise\n3 pTM pTM global\n#\n\
loop_\n_ma_qa_metric_local.model_id\n_ma_qa_metric_local.label_asym_id\n\
_ma_qa_metric_local.label_seq_id\n_ma_qa_metric_local.label_comp_id\n\
_ma_qa_metric_local.metric_id\n_ma_qa_metric_local.metric_value\n\
1 A 1 GLY 1 95.0\n#\n\
loop_\n_ma_qa_metric_local_pairwise.model_id\n\
_ma_qa_metric_local_pairwise.label_asym_id_1\n\
_ma_qa_metric_local_pairwise.label_seq_id_1\n\
_ma_qa_metric_local_pairwise.label_asym_id_2\n\
_ma_qa_metric_local_pairwise.label_seq_id_2\n\
_ma_qa_metric_local_pairwise.metric_id\n\
_ma_qa_metric_local_pairwise.metric_value\n\
1 A 1 A 2 2 3.2\n#\n\
loop_\n_ma_qa_metric_global.model_id\n_ma_qa_metric_global.metric_id\n\
_ma_qa_metric_global.metric_value\n\
1 3 0.82\n#\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n";

const XYZ: &str = "3\nframe-0\nC 0 0 0\nN 1 0 0\nO 0 1 0\n\
3\nframe-1\nC 0 0 0\nN 1 0 0\nO 0 2 0\n";

pub(super) fn register(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    bench_gw_004(group);
    bench_gw_006(group);
    bench_gw_008(group);
    bench_gw_012(group);
    bench_gw_013(group);
    bench_gw_014(group);
    bench_gw_022(group);
    bench_gw_023(group);
    bench_gw_028(group);
    bench_gw_030(group);
    bench_gw_033(group);
    bench_gw_035(group);
    bench_gw_039(group);
    bench_gw_040(group);
}

fn bench_gw_004(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let input = InputBuffer::from_bytes(UNKNOWN_CATEGORY_CIF.as_bytes().to_vec());
    let (document, structure, findings) =
        match molframe::cif::read_with_document(&input, &ReadOptions::new()) {
            Ok(result) => result,
            Err(findings) => panic!("GW-004 setup failed: {findings:?}"),
        };
    assert!(
        findings.is_empty(),
        "GW-004 setup emitted findings: {findings:?}"
    );
    let selection = AtomSelection::All(structure.atom_count());
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-004", |b| {
        b.iter(|| {
            let preserved = molframe::write_preserving(&document);
            let moved = match molframe::transform(
                &structure,
                &selection,
                &Rigid::translation([2.0, 0.0, 0.0]),
            ) {
                Ok(moved) => moved,
                Err(findings) => panic!("GW-004 transform failed: {findings:?}"),
            };
            black_box((preserved.len(), moved.generation().get()));
        });
    });
}

fn bench_gw_006(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = super::fixture(MODEL_CIF, "model.cif");
    group.bench_function("GW-006", |b| {
        b.iter(|| {
            let Some(confidence) = structure.confidence() else {
                panic!("GW-006 confidence metrics are missing");
            };
            black_box((
                confidence.plddt().count(),
                confidence.pae().count(),
                confidence.ptm().count(),
            ));
        });
    });
}

fn bench_gw_008(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    use molframe::QueryStructure as _;

    let structure = structure(Sample::Tiny);
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-008", |b| {
        b.iter(|| {
            black_box(
                structure
                    .select_text(
                        "within 4 of element C",
                        &AnalysisPolicy::default(),
                        &molframe::ExecutionContext::default(),
                    )
                    .is_ok(),
            );
        });
    });
}

fn bench_gw_012(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    group.throughput(Throughput::Elements(structure.residue_count() as u64));
    group.bench_function("GW-012", |b| {
        b.iter(|| {
            let narrow = match molframe::analysis::residue_contact_map(
                &structure,
                4.0,
                1,
                molframe::SpatialBackend::Auto,
                &molframe::ExecutionContext::default(),
            ) {
                Ok(map) => map,
                Err(error) => panic!("GW-012 narrow map failed: {error}"),
            };
            let broad = match molframe::analysis::residue_contact_map(
                &structure,
                8.0,
                1,
                molframe::SpatialBackend::Auto,
                &molframe::ExecutionContext::default(),
            ) {
                Ok(map) => map,
                Err(error) => panic!("GW-012 broad map failed: {error}"),
            };
            black_box((narrow.contacts().len(), broad.contacts().len()));
        });
    });
}

fn bench_gw_013(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    let positions = coordinates(&structure);
    let radii = vec![1.7; positions.len()];
    group.throughput(Throughput::Elements(positions.len() as u64));
    group.bench_function("GW-013", |b| {
        b.iter(|| {
            black_box(
                molframe::surface::shrake_rupley(
                    &positions,
                    &radii,
                    1.4,
                    96,
                    &molframe::ExecutionContext::default(),
                )
                .is_ok(),
            );
        });
    });
}

fn bench_gw_014(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let positions = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]];
    let radii = [1.5, 1.5];
    group.throughput(Throughput::Elements(positions.len() as u64));
    group.bench_function("GW-014", |b| {
        b.iter(|| {
            let buried = molframe::surface::buried_surface(
                &positions,
                &radii,
                1.4,
                96,
                &[true, false],
                &molframe::ExecutionContext::default(),
            );
            let ses = molframe::surface::solvent_excluded_surface(&positions, &radii, 1.4, 1.0);
            black_box((buried.is_ok(), ses.is_ok()));
        });
    });
}

fn bench_gw_022(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let Some(frames) = molframe::traj::parse_xyz(XYZ) else {
        panic!("GW-022 setup XYZ is invalid");
    };
    group.throughput(Throughput::Bytes(XYZ.len() as u64));
    group.bench_function("GW-022", |b| {
        b.iter(|| {
            let written = molframe::traj::write_xyz(&frames);
            let Some(round_trip) = molframe::traj::parse_xyz(&written) else {
                panic!("GW-022 XYZ round trip failed");
            };
            let timesteps: Vec<_> = round_trip
                .iter()
                .enumerate()
                .map(|(frame, value)| molframe::traj::Timestep {
                    frame,
                    positions: value.atoms.iter().map(|atom| atom.position).collect(),
                    ..Default::default()
                })
                .collect();
            let rmsd = match molframe::traj::rmsd_to_reference(
                &timesteps,
                0,
                molframe::traj::FrameAlignment::None,
            ) {
                Ok(rmsd) => rmsd,
                Err(error) => panic!("GW-022 RMSD failed: {error}"),
            };
            black_box((round_trip.len(), rmsd));
        });
    });
}

fn bench_gw_023(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let frames = vec![
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0.1, 0.0, 0.0], [1.1, 0.0, 0.0], [0.0, 1.2, 0.0]],
        vec![[0.2, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 1.4, 0.0]],
    ];
    let views: Vec<_> = frames.iter().map(Vec::as_slice).collect();
    group.throughput(Throughput::Elements(frames.len() as u64));
    group.bench_function("GW-023", |b| {
        b.iter(|| {
            let fluctuation = molframe::geom::rmsf(&views);
            let pca = molframe::traj::cartesian_pca(
                &frames,
                molframe::traj::CartesianFit::None,
                2,
                1024 * 1024,
            );
            black_box((fluctuation.is_ok(), pca.is_ok()));
        });
    });
}

fn bench_gw_028(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    let coordinates = coordinates(&structure);
    group.throughput(Throughput::Elements(coordinates.len() as u64));
    group.bench_function("GW-028", |b| {
        b.iter(|| {
            let lddt = molframe::compare::lddt(
                &coordinates,
                &coordinates,
                15.0,
                &molframe::core::ExecutionContext::default(),
            );
            let tm = molframe::compare::tm_score(&coordinates, &coordinates);
            let ts = molframe::compare::gdt_ts(&coordinates, &coordinates);
            let ha = molframe::compare::gdt_ha(&coordinates, &coordinates);
            black_box((lddt.is_ok(), tm.is_ok(), ts.is_ok(), ha.is_ok()));
        });
    });
}

fn bench_gw_030(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-030", |b| {
        b.iter(|| {
            let first_flags = molframe::validate::quality_flags(&structure);
            let second_flags = molframe::validate::quality_flags(&structure);
            let complete = molframe::validate::completeness(&structure, molframe::Namespace::Label);
            let clashes = molframe::validate::clashes(
                &structure,
                0.4,
                molframe::RadiusSet::Bondi,
                molframe::SpatialBackend::Auto,
                &molframe::ExecutionContext::default(),
            );
            black_box((
                first_flags.len(),
                second_flags.len(),
                complete.is_ok(),
                clashes.is_ok(),
            ));
        });
    });
}

fn bench_gw_033(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    let table = molframe::AtomArrowTable::new(&structure);
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-033/arrow_stream", |b| {
        b.iter(|| black_box(table.arrow_stream()));
    });
}

fn bench_gw_035(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let structure = structure(Sample::Tiny);
    let options = molframe::GraphOptions {
        nodes: molframe::NodeLevel::Atoms,
        edges: molframe::EdgeKind::Radius { cutoff: 3.0 },
        direction: molframe::EdgeDirection::Symmetric,
        node_features: vec![
            molframe::NodeFeature::PositionX,
            molframe::NodeFeature::PositionY,
            molframe::NodeFeature::PositionZ,
        ],
        edge_features: vec![molframe::EdgeFeature::Distance],
        missing: molframe::MissingFeaturePolicy::Error,
        backend: molframe::SpatialBackend::Auto,
        periodic: false,
    };
    group.throughput(Throughput::Elements(structure.atom_count().into()));
    group.bench_function("GW-035", |b| {
        b.iter(|| {
            black_box(molframe::graph(
                &structure,
                &options,
                &molframe::ExecutionContext::default(),
            ))
        });
    });
}

fn bench_gw_039(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let space = molframe::PolicySpace::new(AnalysisPolicy::default())
        .vary(molframe::PolicyDimension::altloc([
            AltlocPolicy::KeepAll,
            AltlocPolicy::First,
        ]))
        .vary(molframe::PolicyDimension::model([
            ModelChoice::First,
            ModelChoice::All,
        ]));
    group.bench_function("GW-039", |b| {
        b.iter(|| black_box(space.clone().plan()));
    });
}

fn bench_gw_040(group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>) {
    let profile = molframe::fx::motifbench_1_0();
    let metrics = BTreeMap::from([("rmsd".into(), 1.0), ("motif_rmsd".into(), 0.5)]);
    group.bench_function("GW-040", |b| {
        b.iter(|| black_box(profile.decide_candidate(&metrics)));
    });
}
