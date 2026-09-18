#[cfg(feature = "compare")]
use super::{ComparisonMetric, ComparisonRequest};
use super::{
    ContactsRequest, CoordinateInput, CoordinateSlot, ExecutionPlanError, FrameInput,
    GeometryRequest, Plan, PlanOperation, PlanValue, RmsdRequest, ScalarInput, SelectionRequest,
    SpatialRequest, StructureRequest, StructureValue,
};

mod physical_tests;
mod structure_tests;
#[cfg(feature = "surface")]
mod surface_tests;
#[cfg(feature = "traj")]
mod trajectory_tests;
use crate::BondInference;
use crate::{
    AnalysisPolicy, ReadOptions, SpatialBackend, SpatialSearchOptions, Status, read_bytes,
};

const STRUCTURE_CIF: &str = "data_operation\n\
loop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_entity_id\n\
_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n\
_atom_site.auth_seq_id\n\
_atom_site.auth_comp_id\n\
_atom_site.auth_asym_id\n\
_atom_site.auth_atom_id\n\
_atom_site.pdbx_PDB_model_num\n\
ATOM 1 C CA . GLY A 1 1 ? 0 0 0 1.0 10.0 1 GLY A CA 1\n";

#[test]
fn operation_requests_validate_before_plan_execution() {
    let invalid = ContactsRequest::new(
        "",
        "B",
        4.5,
        SpatialBackend::Auto,
        AnalysisPolicy::default(),
    );
    assert!(matches!(
        invalid,
        Err(ExecutionPlanError::InvalidRequest(_))
    ));

    let request = ContactsRequest::new(
        "chain A",
        "chain B",
        4.5,
        SpatialBackend::Auto,
        AnalysisPolicy::default(),
    )
    .expect("valid request");
    let mut plan = Plan::new();
    plan.add("contacts", PlanOperation::Contacts(Box::new(request)))
        .expect("first operation");
    let duplicate = plan.add(
        "contacts",
        RmsdRequest::new(CoordinateSlot::new(0), CoordinateSlot::new(1)),
    );
    assert!(matches!(duplicate, Err(ExecutionPlanError::DuplicateId(_))));
}

#[test]
fn array_operations_execute_in_native_sorted_order() {
    let mobile = [[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let reference = mobile;
    let arrays = [
        super::CoordinateInput { positions: &mobile },
        super::CoordinateInput {
            positions: &reference,
        },
    ];
    let mut plan = Plan::new();
    plan.add(
        "zeta",
        RmsdRequest::new(CoordinateSlot::new(0), CoordinateSlot::new(1)),
    )
    .expect("first RMSD");
    plan.add(
        "alpha",
        RmsdRequest::new(CoordinateSlot::new(0), CoordinateSlot::new(1)),
    )
    .expect("second RMSD");

    let result = plan
        .execute(
            super::PlanInput {
                structure: None,
                arrays: &arrays,
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native plan");
    let ids = result
        .entries
        .iter()
        .map(|entry| entry.id.as_ref())
        .collect::<Vec<_>>();
    assert_eq!(ids, ["alpha", "zeta"]);
    assert!(result.entries.iter().all(|entry| matches!(
        entry.value,
        super::PlanValue::Rmsd(value) if value == 0.0
    )));
}

#[test]
fn geometry_operations_match_direct_kernels_without_python_data_loops() {
    let positions = [
        [-1.0_f32, 0.0, 0.0],
        [1.0_f32, 0.0, 0.0],
        [0.0_f32, 1.0, 0.0],
    ];
    let masses = [1.0_f64, 2.0, 3.0];
    let frame_a = positions;
    let frame_b = [[-0.5_f32, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]];
    let frame_positions = [
        frame_a[0], frame_a[1], frame_a[2], frame_b[0], frame_b[1], frame_b[2],
    ];
    let arrays = [CoordinateInput {
        positions: &positions,
    }];
    let scalars = [ScalarInput { values: &masses }];
    let frames = [FrameInput {
        positions: &frame_positions,
        frame_count: 2,
        atom_count: positions.len(),
    }];

    let mut plan = Plan::new();
    plan.add(
        "centroid",
        PlanOperation::Geometry(Box::new(GeometryRequest::Centroid { positions: 0 })),
    )
    .expect("centroid operation");
    plan.add(
        "radius",
        PlanOperation::Geometry(Box::new(GeometryRequest::RadiusOfGyration {
            positions: 0,
            masses: Some(0),
        })),
    )
    .expect("radius operation");
    plan.add(
        "matrix",
        PlanOperation::Geometry(Box::new(GeometryRequest::DistanceMatrix { positions: 0 })),
    )
    .expect("distance matrix operation");
    plan.add(
        "rmsf",
        PlanOperation::Geometry(Box::new(GeometryRequest::Rmsf { frames: 0 })),
    )
    .expect("RMSF operation");

    let result = plan
        .execute(
            super::PlanInput {
                structure: None,
                arrays: &arrays,
                scalars: &scalars,
                frames: &frames,
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native geometry plan");

    let direct_centroid = molframe_geom::centroid(&positions);
    let direct_radius = molframe_geom::radius_of_gyration(&positions, &masses);
    let direct_matrix = molframe_geom::distance_matrix(&positions).expect("direct matrix");
    let direct_rmsf = molframe_geom::rmsf(&[&frame_a, &frame_b]).expect("direct RMSF");

    let values = result
        .entries
        .into_iter()
        .map(|entry| (entry.id.to_string(), entry.value))
        .collect::<std::collections::BTreeMap<_, _>>();
    assert!(matches!(
        values.get("centroid"),
        Some(PlanValue::Geometry(value))
            if matches!(value.as_ref(), super::GeometryValue::Centroid(actual) if *actual == direct_centroid)
    ));
    assert!(matches!(
        values.get("radius"),
        Some(PlanValue::Geometry(value))
            if matches!(value.as_ref(), super::GeometryValue::RadiusOfGyration(actual) if *actual == direct_radius)
    ));
    assert!(matches!(
        values.get("matrix"),
        Some(PlanValue::Geometry(value))
            if matches!(value.as_ref(), super::GeometryValue::DistanceMatrix(actual) if actual == &direct_matrix)
    ));
    assert!(matches!(
        values.get("rmsf"),
        Some(PlanValue::Geometry(value))
            if matches!(value.as_ref(), super::GeometryValue::Rmsf(actual) if actual == &direct_rmsf)
    ));
}

#[test]
fn selection_operations_match_the_facade_query_kernel() {
    let (structure, _) = read_bytes(
        STRUCTURE_CIF.as_bytes().to_vec(),
        Some("operation.cif"),
        &ReadOptions::new(),
    )
    .expect("structure fixture");
    let request = SelectionRequest::new("name CA", AnalysisPolicy::default())
        .expect("compiled selection request");
    let direct =
        crate::QueryStructure::select_text(&structure, "name CA", &AnalysisPolicy::default())
            .expect("direct query selection");
    let mut plan = Plan::new();
    plan.add("selected", request).expect("selection operation");
    let result = plan
        .execute(
            super::PlanInput {
                structure: Some(&structure),
                arrays: &[],
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native selection plan");
    let Some(PlanValue::Selection(value)) = result.entries.first().map(|entry| &entry.value) else {
        panic!("unexpected selection result type");
    };
    assert_eq!(value.selection, direct.selection);
    assert_eq!(value.warnings, direct.warnings);
}

#[test]
fn spatial_operations_match_the_facade_spatial_kernel() {
    let positions = [
        [0.0_f32, 0.0, 0.0],
        [1.0_f32, 0.0, 0.0],
        [0.0_f32, 3.0, 0.0],
    ];
    let arrays = [CoordinateInput {
        positions: &positions,
    }];
    let left = crate::AtomSelection::from_sorted(vec![0, 1]);
    let right = crate::AtomSelection::from_sorted(vec![1, 2]);
    let options = SpatialSearchOptions::with_backend(SpatialBackend::CellList);
    let (direct_pairs, direct_again, direct_within) =
        direct_spatial_results(&positions, &left, &right, options);

    let mut plan = Plan::new();
    plan.add(
        "pairs",
        PlanOperation::Spatial(Box::new(SpatialRequest::NeighborPairs {
            positions: 0,
            left,
            right: right.clone(),
            cutoff: 1.5,
            options,
            periodic: None,
        })),
    )
    .expect("pairs operation");
    plan.add(
        "pairs_again",
        PlanOperation::Spatial(Box::new(SpatialRequest::NeighborPairs {
            positions: 0,
            left: crate::AtomSelection::from_sorted(vec![0]),
            right: right.clone(),
            cutoff: 1.5,
            options,
            periodic: None,
        })),
    )
    .expect("repeated pairs operation");
    plan.add(
        "within",
        PlanOperation::Spatial(Box::new(SpatialRequest::AtomsWithin {
            positions: 0,
            query: crate::AtomSelection::All(3),
            target: crate::AtomSelection::from_sorted(vec![0]),
            cutoff: 1.5,
            options,
            periodic: None,
        })),
    )
    .expect("within operation");

    let result = plan
        .execute(
            super::PlanInput {
                arrays: &arrays,
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("spatial plan");
    assert_eq!(result.cached_index_count, 2);
    let values = result
        .entries
        .into_iter()
        .map(|entry| (entry.id.to_string(), entry.value))
        .collect::<std::collections::BTreeMap<_, _>>();
    assert!(matches!(
        values.get("pairs"),
        Some(PlanValue::Spatial(value))
            if matches!(value.as_ref(), super::SpatialValue::NeighborPairs(actual) if actual == &direct_pairs)
    ));
    assert!(matches!(
        values.get("within"),
        Some(PlanValue::Spatial(value))
            if matches!(value.as_ref(), super::SpatialValue::AtomsWithin(actual) if actual == &direct_within)
    ));
    assert!(matches!(
        values.get("pairs_again"),
        Some(PlanValue::Spatial(value))
            if matches!(value.as_ref(), super::SpatialValue::NeighborPairs(actual) if actual == &direct_again)
    ));
}

fn direct_spatial_results(
    positions: &[[f32; 3]],
    left: &crate::AtomSelection,
    right: &crate::AtomSelection,
    options: SpatialSearchOptions,
) -> (
    Vec<crate::NeighborPair>,
    Vec<crate::NeighborPair>,
    crate::AtomSelection,
) {
    let context = molframe_core::ExecutionContext::default();
    let pairs =
        crate::pairs_within_with_options(positions, left, right, 1.5, options, None, &context)
            .expect("direct spatial pairs");
    let origin = crate::AtomSelection::from_sorted(vec![0]);
    let again =
        crate::pairs_within_with_options(positions, &origin, right, 1.5, options, None, &context)
            .expect("direct repeated spatial pairs");
    let within = crate::within_with_options(
        positions,
        &crate::AtomSelection::All(3),
        &origin,
        1.5,
        options,
        None,
        &context,
    )
    .expect("direct spatial within");
    (pairs, again, within)
}

#[test]
fn chemistry_operations_use_the_same_native_kernel_as_direct_calls() {
    let (structure, _) = read_bytes(
        STRUCTURE_CIF.as_bytes().to_vec(),
        Some("operation.cif"),
        &ReadOptions::new(),
    )
    .expect("structure fixture");
    let options = BondInference::default();
    let direct = super::super::infer_bonds(
        &structure,
        options,
        &molframe_core::ExecutionContext::default(),
    )
    .expect("direct inference");
    let mut plan = Plan::new();
    plan.add("bonds", options)
        .expect("bond inference operation");
    let result = plan
        .execute(
            super::PlanInput {
                structure: Some(&structure),
                arrays: &[],
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native chemistry plan");
    let Some(entry) = result.entries.first() else {
        panic!("missing bond inference result");
    };
    let PlanValue::BondInference(report) = &entry.value else {
        panic!("unexpected bond inference result type");
    };
    assert_eq!(report.structure.atom_count(), direct.structure.atom_count());
    assert_eq!(report.skipped_atoms, direct.skipped_atoms);
    assert_eq!(
        report.structure.data().bonds.iter().count(),
        direct.structure.data().bonds.iter().count()
    );
}

#[cfg(feature = "validate")]
#[test]
fn structure_operations_return_governed_native_analysis() {
    let (structure, _) = read_bytes(
        STRUCTURE_CIF.as_bytes().to_vec(),
        Some("operation.cif"),
        &ReadOptions::new(),
    )
    .expect("structure fixture");
    let mut plan = Plan::new();
    plan.add(
        "quality",
        PlanOperation::Structure(Box::new(StructureRequest::Quality {
            policy: AnalysisPolicy::default(),
        })),
    )
    .expect("quality operation");

    let result = plan
        .execute(
            super::PlanInput {
                structure: Some(&structure),
                arrays: &[],
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native structure operation");
    let Some(entry) = result.entries.first() else {
        panic!("missing structure result");
    };
    let PlanValue::Structure(value) = &entry.value else {
        panic!("unexpected result type");
    };
    let StructureValue::Quality(analysis) = value.as_ref() else {
        panic!("unexpected structure result type");
    };
    assert_eq!(analysis.status, Status::Complete);
    assert!(analysis.value.is_empty());
    assert_eq!(analysis.coverage.intended, structure.atom_count());
}

#[cfg(feature = "compare")]
#[test]
fn comparison_operations_share_borrowed_coordinate_inputs() {
    let mobile = [[0.0_f32, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let reference = mobile;
    let arrays = [
        super::CoordinateInput { positions: &mobile },
        super::CoordinateInput {
            positions: &reference,
        },
    ];
    let mut plan = Plan::new();
    plan.add(
        "lddt",
        ComparisonRequest::new(
            CoordinateSlot::new(0),
            CoordinateSlot::new(1),
            ComparisonMetric::Lddt {
                inclusion_radius: 15.0,
            },
        ),
    )
    .expect("lDDT operation");
    plan.add(
        "tm",
        ComparisonRequest::new(
            CoordinateSlot::new(0),
            CoordinateSlot::new(1),
            ComparisonMetric::TmScore,
        ),
    )
    .expect("TM-score operation");
    plan.add(
        "gdt",
        ComparisonRequest::new(
            CoordinateSlot::new(0),
            CoordinateSlot::new(1),
            ComparisonMetric::GdtTs,
        ),
    )
    .expect("GDT operation");

    let result = plan
        .execute(
            super::PlanInput {
                structure: None,
                arrays: &arrays,
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("native comparison plan");
    assert_eq!(result.entries.len(), 3);
    assert!(result.entries.iter().all(|entry| {
        matches!(
            entry.value,
            super::PlanValue::Comparison(value) if (value.value - 1.0).abs() < f64::EPSILON
        )
    }));
}
