//! Executable golden fixtures for the native analysis and interoperability rows.

use std::collections::BTreeMap;

use molframe::{
    AnalysisPolicy, AnnotationColumn, AtomAnnotation, AtomIndex, AtomSelection, BondOrder,
    BondProvenance, BondRecord, BondTableBuilder, InputBuffer, ModelCifExt, Presence, ReadOptions,
    Rigid,
};
use molframe_bench::{Sample, coordinates, structure};

#[path = "golden/native_pending.rs"]
mod native_pending;

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
loop_\n_ma_qa_metric_global.model_id\n_ma_qa_metric_global.metric_id\n_ma_qa_metric_global.metric_value\n\
1 3 0.82\n#\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C CA GLY A 1 0 0 0\n";

const RAMA_CIF: &str = "data_rama\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N ALA A 1 0 0 0\n\
ATOM 2 C CA ALA A 1 1.46 0 0\n\
ATOM 3 C C ALA A 1 2.0 1.2 0\n\
ATOM 4 N N ALA A 2 3.3 1.4 0\n\
ATOM 5 C CA ALA A 2 4.0 2.5 0\n\
ATOM 6 C C ALA A 2 5.4 2.5 0\n\
ATOM 7 N N ALA A 3 6.0 3.6 0\n\
ATOM 8 C CA ALA A 3 7.0 3.6 0\n\
ATOM 9 C C ALA A 3 8.0 4.0 0\n";

fn read_fixture(text: &str, name: &str) -> molframe::Structure {
    match molframe::read_bytes(text.as_bytes().to_vec(), Some(name), &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("golden read failed: {findings:?}"),
    }
}

#[test]
fn gw_004_preserves_unknown_categories_while_coordinates_are_edited() {
    let input = InputBuffer::from_bytes(UNKNOWN_CATEGORY_CIF.as_bytes().to_vec());
    let (document, structure, findings) =
        match molframe::cif::read_with_document(&input, &ReadOptions::new()) {
            Ok(result) => result,
            Err(findings) => panic!("document read failed: {findings:?}"),
        };
    assert!(findings.is_empty());
    assert!(molframe::write_preserving(&document).contains("_custom.note"));

    let moved = molframe::transform(
        &structure,
        &AtomSelection::All(structure.atom_count()),
        &Rigid::translation([2.0, 0.0, 0.0]),
    )
    .unwrap_or_else(|findings| panic!("coordinate edit failed: {findings:?}"));
    assert!(
        moved.positions()[0]
            .iter()
            .zip([2.0_f32, 0.0, 0.0])
            .all(|(actual, expected)| (*actual - expected).abs() < 1.0e-6)
    );
}

#[test]
fn gw_006_extracts_modelcif_plddt_pae_and_ptm_metrics() {
    let structure = read_fixture(MODEL_CIF, "model.cif");
    let confidence = structure
        .confidence()
        .unwrap_or_else(|| panic!("ModelCIF confidence is missing"));
    assert_eq!(confidence.plddt().count(), 1);
    assert_eq!(confidence.pae().count(), 1);
    assert_eq!(confidence.ptm().count(), 1);
    assert_eq!(
        confidence.plddt().next().map(|metric| metric.value),
        Some(95.0)
    );
}

#[test]
fn gw_008_evaluates_a_spatial_binding_pocket_selection() {
    use molframe::QueryStructure as _;

    let structure = structure(Sample::Tiny);
    let evaluation = structure
        .select_text(
            "within 4 of element C",
            &AnalysisPolicy::default(),
            &molframe::ExecutionContext::default(),
        )
        .unwrap_or_else(|findings| panic!("selection failed: {findings:?}"));
    let selected = evaluation.selection.iter().count();
    assert!(selected > 0);
    assert!(evaluation.warnings.is_empty());
}

#[test]
fn gw_010_classifies_backbone_torsions_against_versioned_reference_data() {
    let structure = rama_structure();
    let grid = molframe::validate::ReferenceDistribution::grid(
        "general",
        vec![-180.0, 0.0, 180.0],
        vec![-180.0, 0.0, 180.0],
        vec![1.0, 1.0, 1.0, 1.0],
    )
    .unwrap_or_else(|error| panic!("Ramachandran grid failed: {error}"));
    let library = molframe::validate::ReferenceLibrary::new("rama", "golden-1", [grid])
        .unwrap_or_else(|error| panic!("Ramachandran library failed: {error}"));
    let basin = molframe::validate::RamachandranBasin::new(
        molframe::validate::RamachandranRegion::AlphaHelixRight,
        "general",
    )
    .unwrap_or_else(|error| panic!("Ramachandran basin failed: {error}"));
    let options = molframe::validate::RamachandranOptions::new(&library, [basin], 0.0)
        .unwrap_or_else(|error| panic!("Ramachandran options failed: {error}"));
    let records = molframe::validate::ramachandran(&structure, &options)
        .unwrap_or_else(|error| panic!("Ramachandran workflow failed: {error}"));
    assert_eq!(records.len(), 1);
    assert!(records.iter().all(|record| {
        record.phi.is_finite()
            && record.psi.is_finite()
            && record.region == molframe::validate::RamachandranRegion::AlphaHelixRight
            && record.reference.version.as_ref() == "golden-1"
    }));
    assert!(
        molframe::validate::ramachandran_outliers(&structure, &options)
            .unwrap_or_else(|error| panic!("Ramachandran outlier workflow failed: {error}"))
            .is_empty()
    );
}

#[test]
fn gw_012_contact_maps_are_monotonic_across_cutoffs() {
    let structure = structure(Sample::Tiny);
    let narrow = molframe::analysis::residue_contact_map(
        &structure,
        4.0,
        1,
        molframe::SpatialBackend::Auto,
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("narrow contact map failed: {error}"));
    let broad = molframe::analysis::residue_contact_map(
        &structure,
        8.0,
        1,
        molframe::SpatialBackend::Auto,
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("broad contact map failed: {error}"));
    assert_eq!(narrow.residue_count(), structure.residue_count());
    assert!(broad.contacts().len() >= narrow.contacts().len());
    assert!(
        broad
            .contacts()
            .windows(2)
            .all(|pair| (pair[0].first, pair[0].second) <= (pair[1].first, pair[1].second))
    );
}

#[test]
fn gw_013_reports_per_atom_sasa_with_finite_areas() {
    let structure = structure(Sample::Tiny);
    let positions = coordinates(&structure);
    let radii = vec![1.7; positions.len()];
    let areas = molframe::surface::shrake_rupley(
        &positions,
        &radii,
        1.4,
        96,
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("SASA failed: {error}"));
    assert_eq!(areas.len(), positions.len());
    assert!(areas.iter().all(|area| area.is_finite() && *area >= 0.0));
    assert!(areas.iter().sum::<f64>() > 0.0);
}

#[test]
fn gw_014_agrees_on_buried_and_solvent_excluded_surface() {
    let positions = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]];
    let radii = [1.5, 1.5];
    let buried = molframe::surface::buried_surface(
        &positions,
        &radii,
        1.4,
        96,
        &[true, false],
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("buried surface failed: {error}"));
    assert!(buried.buried > 0.0);
    assert!(buried.first_alone + buried.second_alone >= buried.together);

    let ses = molframe::surface::solvent_excluded_surface(&positions, &radii, 1.4, 1.0)
        .unwrap_or_else(|error| panic!("SES failed: {error}"));
    assert!(!ses.triangles.is_empty());
    assert!(ses.area.is_finite() && ses.area > 0.0);
}

#[test]
fn gw_022_round_trips_xyz_and_measures_streamed_rmsd() {
    let source = "3\nframe-0\nC 0 0 0\nN 1 0 0\nO 0 1 0\n\
3\nframe-1\nC 0 0 0\nN 1 0 0\nO 0 2 0\n";
    let frames = molframe::traj::parse_xyz(source).unwrap_or_else(|| panic!("XYZ parse failed"));
    let written = molframe::traj::write_xyz(&frames);
    let round_trip =
        molframe::traj::parse_xyz(&written).unwrap_or_else(|| panic!("XYZ reparse failed"));
    assert_eq!(round_trip, frames);

    let timesteps: Vec<_> = frames
        .iter()
        .enumerate()
        .map(|(frame, value)| molframe::traj::Timestep {
            frame,
            positions: value.atoms.iter().map(|atom| atom.position).collect(),
            ..Default::default()
        })
        .collect();
    let rmsd =
        molframe::traj::rmsd_to_reference(&timesteps, 0, molframe::traj::FrameAlignment::None)
            .unwrap_or_else(|error| panic!("trajectory RMSD failed: {error}"));
    assert_eq!(rmsd.len(), 2);
    assert!(rmsd[0].abs() < 1.0e-6);
    assert!(rmsd[1] > 0.0);

    let directory = tempfile::tempdir().expect("trajectory directory");
    let path = directory.path().join("frames.trr");
    let encoded = molframe::traj::write_trr(&timesteps, molframe::traj::TrrWriteOptions::default())
        .expect("TRR encode");
    std::fs::write(&path, encoded).expect("TRR fixture");
    let mut reader =
        molframe::traj::read_trajectory(&path, &molframe::traj::TrajectoryReaderOptions::default())
            .expect("pull reader");
    let context = molframe::ExecutionContext::default();
    let mut streamed = Vec::new();
    molframe::traj::rmsd_stream(
        &mut *reader,
        &timesteps[0].positions,
        molframe::traj::FrameAlignment::None,
        &context,
        16384,
        |_, _, value| {
            streamed.push(value);
            Ok(())
        },
    )
    .expect("streamed RMSD");
    assert_eq!(streamed, rmsd);
    assert_eq!(context.reserved_bytes(), 0);
}

#[test]
fn gw_023_rmsf_and_cartesian_pca_are_deterministic() {
    let frames = vec![
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0.1, 0.0, 0.0], [1.1, 0.0, 0.0], [0.0, 1.2, 0.0]],
        vec![[0.2, 0.0, 0.0], [1.2, 0.0, 0.0], [0.0, 1.4, 0.0]],
    ];
    let views: Vec<_> = frames.iter().map(Vec::as_slice).collect();
    let fluctuation =
        molframe::geom::rmsf(&views).unwrap_or_else(|error| panic!("RMSF failed: {error:?}"));
    assert_eq!(fluctuation.len(), 3);
    assert!(fluctuation.iter().any(|value| *value > 0.0));

    let pca =
        molframe::traj::cartesian_pca(&frames, molframe::traj::CartesianFit::None, 2, 1024 * 1024)
            .unwrap_or_else(|error| panic!("Cartesian PCA failed: {error}"));
    assert_eq!(pca.eigenvalues.len(), 2);
    assert_eq!(pca.projections.len(), frames.len());
}

#[test]
fn gw_028_structure_scores_are_perfect_for_identical_coordinates() {
    let structure = structure(Sample::Tiny);
    let coordinates = coordinates(&structure);
    let lddt = molframe::compare::lddt(
        &coordinates,
        &coordinates,
        15.0,
        &molframe::core::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("lDDT failed: {error}"));
    let tm = molframe::compare::tm_score(&coordinates, &coordinates)
        .unwrap_or_else(|error| panic!("TM-score failed: {error}"));
    let ts = molframe::compare::gdt_ts(&coordinates, &coordinates)
        .unwrap_or_else(|error| panic!("GDT-TS failed: {error}"));
    let ha = molframe::compare::gdt_ha(&coordinates, &coordinates)
        .unwrap_or_else(|error| panic!("GDT-HA failed: {error}"));
    for score in [lddt, tm, ts, ha] {
        assert!((score - 1.0).abs() < 1.0e-6);
    }
}

#[test]
fn gw_030_validation_report_is_stable_and_structured() {
    let structure = structure(Sample::Tiny);
    let first_flags = molframe::validate::quality_flags(&structure);
    let second_flags = molframe::validate::quality_flags(&structure);
    assert_eq!(first_flags, second_flags);
    assert!(molframe::validate::completeness(&structure, molframe::Namespace::Label).is_ok());
    let clashes = molframe::validate::clashes(
        &structure,
        0.4,
        molframe::RadiusSet::Bondi,
        molframe::SpatialBackend::Auto,
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("clash report failed: {error}"));
    assert!(clashes.iter().all(|clash| clash.overlap.is_finite()));
}

#[test]
fn gw_033_writes_a_nonempty_arrow_atom_table() {
    let structure = structure(Sample::Tiny);
    let file = tempfile::NamedTempFile::new()
        .unwrap_or_else(|error| panic!("temporary Arrow file failed: {error}"));
    molframe::write_atom_ipc(file.path(), &structure)
        .unwrap_or_else(|error| panic!("Arrow IPC export failed: {error}"));
    let size = file
        .as_file()
        .metadata()
        .unwrap_or_else(|error| panic!("Arrow metadata failed: {error}"))
        .len();
    assert!(size > 0);
}

#[test]
fn gw_035_radius_graph_has_typed_nodes_and_edges() {
    let structure = structure(Sample::Tiny);
    let graph = molframe::graph(
        &structure,
        &molframe::GraphOptions {
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
        },
        &molframe::ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("radius graph failed: {error}"));
    assert_eq!(graph.node_count, structure.atom_count() as usize);
    assert_eq!(graph.node_features.columns, 3);
    assert_eq!(graph.edge_features.columns, 1);
    assert!(graph.edge_count > 0);
}

#[test]
fn gw_039_policy_audit_expands_in_deterministic_order() {
    let space = molframe::PolicySpace::new(AnalysisPolicy::default())
        .vary(molframe::PolicyDimension::altloc([
            molframe::AltlocPolicy::KeepAll,
            molframe::AltlocPolicy::First,
        ]))
        .vary(molframe::PolicyDimension::model([
            molframe::ModelChoice::First,
            molframe::ModelChoice::All,
        ]));
    let plan = space
        .plan()
        .unwrap_or_else(|error| panic!("policy plan failed: {error}"));
    assert_eq!(plan.cost(), 4);
    assert_eq!(plan.policies().len(), plan.coordinates().len());
    assert_eq!(plan.fields().len(), 2);
}

#[test]
fn gw_040_versioned_fx_profile_returns_a_stable_verdict() {
    let profile = molframe::fx::motifbench_1_0();
    let metrics = BTreeMap::from([("rmsd".into(), 1.0), ("motif_rmsd".into(), 0.5)]);
    let verdict = profile.decide_candidate(&metrics);
    assert_eq!(profile.id(), "motifbench-1.0");
    assert_eq!(verdict.status, molframe::fx::VerdictStatus::Pass);
}

fn rama_structure() -> molframe::Structure {
    let source = read_fixture(RAMA_CIF, "rama.cif");
    let mut data = source.data().clone();
    let roles = [
        molframe::PolymerAtomRole::PROTEIN_NITROGEN,
        molframe::PolymerAtomRole::PROTEIN_ALPHA_CARBON,
        molframe::PolymerAtomRole::PROTEIN_CARBONYL_CARBON,
    ];
    let role_values = AnnotationColumn::from_entries(
        (0..9).map(|index| (roles[index % roles.len()].code(), Presence::Present)),
    )
    .unwrap_or_else(|error| panic!("Ramachandran role annotation failed: {error}"));
    data.annotations.insert(
        molframe::POLYMER_ATOM_ROLE_ANNOTATION,
        AtomAnnotation::Integer(role_values),
    );
    let mut bonds = BondTableBuilder::new();
    for (carbon, nitrogen) in [(2, 3), (5, 6)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(carbon),
            atom_b: AtomIndex::new(nitrogen),
            order: BondOrder::Single,
            provenance: BondProvenance::User,
        });
    }
    data.bonds = bonds.finish();
    molframe::Structure::new(data)
}
