use crate::{
    AnalysisPolicy, FloatInput, PhysicalRequest, PhysicalValue, Plan, PlanInput, PlanOperation,
    PlanValue, ReadOptions, ScalarInput, SpatialBackend,
};

const PHYSICAL_CIF: &str = r"data_physical
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . GLY A 1 1 ? 0 0 0 1.0 10.0 1 GLY A CA 1
ATOM 2 C CB . GLY A 1 1 ? 1 0 0 1.0 10.0 1 GLY A CB 1
ATOM 3 C CG . GLY A 1 1 ? 2 0 0 1.0 10.0 1 GLY A CG 1
ATOM 4 C CD . GLY A 1 1 ? 20 0 0 1.0 10.0 1 GLY A CD 1
";

fn structure() -> crate::Structure {
    crate::read_bytes(
        PHYSICAL_CIF.as_bytes().to_vec(),
        Some("physical.cif"),
        &ReadOptions::new(),
    )
    .expect("physical fixture")
    .0
}

struct DirectPhysicalAnalyses {
    radial: crate::Analysis<Vec<molframe_analysis::RadialBin>>,
    coordination: crate::Analysis<Vec<u32>>,
    leaflets: crate::Analysis<Vec<molframe_analysis::Leaflet>>,
    linear: crate::Analysis<Vec<molframe_analysis::LinearDensityBin>>,
    density: crate::Analysis<molframe_analysis::DensityGrid>,
    pore: crate::Analysis<Vec<molframe_analysis::PoreSample>>,
    surface_contacts: crate::Analysis<Vec<molframe_analysis::Contact>>,
}

struct PhysicalOptions {
    radial: molframe_analysis::RadialDistributionOptions,
    leaflets: molframe_analysis::LeafletOptions,
    linear: molframe_analysis::LinearDensityOptions,
    grid: molframe_analysis::DensityGridSpec,
    pore: molframe_analysis::PoreProfileOptions,
}

fn options() -> PhysicalOptions {
    PhysicalOptions {
        radial: molframe_analysis::RadialDistributionOptions {
            minimum_distance: 0.0,
            maximum_distance: 3.5,
            bins: 7,
            volume: 1_000.0,
            backend: SpatialBackend::BruteForce,
        },
        leaflets: molframe_analysis::LeafletOptions {
            connection_distance: 1.5,
            backend: SpatialBackend::BruteForce,
        },
        linear: molframe_analysis::LinearDensityOptions {
            axis: molframe_analysis::CartesianAxis::X,
            minimum: -1.0,
            maximum: 22.0,
            bins: 5,
        },
        grid: molframe_analysis::DensityGridSpec {
            origin: [-1.0, -1.0, -1.0],
            spacing: [5.0, 5.0, 5.0],
            shape: [5, 2, 2],
        },
        pore: molframe_analysis::PoreProfileOptions {
            axis_origin: [0.0, 0.0, 0.0],
            axis_direction: [1.0, 0.0, 0.0],
            start: 0.0,
            end: 20.0,
            samples: 4,
            search_radius: 4.0,
            grid_spacing: 1.0,
            probe_radius: 1.0,
            memory_limit_bytes: 100_000_000,
        },
    }
}

fn direct_analyses(
    structure: &crate::Structure,
    left: &crate::AtomSelection,
    weights: &[f64],
    radii: &[f32],
    options: &PhysicalOptions,
    policy: &AnalysisPolicy,
) -> DirectPhysicalAnalyses {
    let radial_kernel = molframe_analysis::radial_distribution_kernel(left, left, options.radial);
    let coordination_kernel = molframe_analysis::coordination_numbers_kernel(
        left,
        left,
        0.0,
        2.1,
        SpatialBackend::BruteForce,
    );
    let leaflets_kernel = molframe_analysis::leaflets_kernel(left, options.leaflets);
    let linear_kernel = molframe_analysis::linear_density_kernel(weights, options.linear);
    let density_kernel = molframe_analysis::density_map_kernel(weights, options.grid);
    let pore_kernel = molframe_analysis::pore_profile_kernel(radii, options.pore);
    let contacts_kernel = molframe_analysis::surface_contacts_kernel(
        radii,
        molframe_analysis::SurfaceContactOptions {
            tolerance: 0.5,
            probe: 1.0,
            surface_density: 32.0,
            minimum_area: 0.2,
            backend: SpatialBackend::BruteForce,
        },
    );
    DirectPhysicalAnalyses {
        radial: molframe_analysis::analyse_structure(
            structure,
            policy,
            &radial_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct radial analysis"),
        coordination: molframe_analysis::analyse_structure(
            structure,
            policy,
            &coordination_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct coordination analysis"),
        leaflets: molframe_analysis::analyse_structure(
            structure,
            policy,
            &leaflets_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct leaflet analysis"),
        linear: molframe_analysis::analyse_structure(
            structure,
            policy,
            &linear_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct linear density analysis"),
        density: molframe_analysis::analyse_structure(
            structure,
            policy,
            &density_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct density analysis"),
        pore: molframe_analysis::analyse_structure(
            structure,
            policy,
            &pore_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct pore analysis"),
        surface_contacts: molframe_analysis::analyse_structure(
            structure,
            policy,
            &contacts_kernel,
            &molframe_core::ExecutionContext::default(),
        )
        .expect("direct surface contacts analysis"),
    }
}

fn add_structure_operations(
    plan: &mut Plan,
    left: &crate::AtomSelection,
    options: &PhysicalOptions,
    policy: &AnalysisPolicy,
) {
    plan.add(
        "radial",
        PlanOperation::Physical(Box::new(PhysicalRequest::RadialDistribution {
            left: left.clone(),
            right: left.clone(),
            options: options.radial,
            periodic: false,
            policy: policy.clone(),
        })),
    )
    .expect("radial operation");
    plan.add(
        "coordination",
        PlanOperation::Physical(Box::new(PhysicalRequest::CoordinationNumbers {
            left: left.clone(),
            right: left.clone(),
            minimum_distance: 0.0,
            maximum_distance: 2.1,
            backend: SpatialBackend::BruteForce,
            periodic: false,
            policy: policy.clone(),
        })),
    )
    .expect("coordination operation");
    plan.add(
        "leaflets",
        PlanOperation::Physical(Box::new(PhysicalRequest::Leaflets {
            sites: left.clone(),
            options: options.leaflets,
            periodic: false,
            policy: policy.clone(),
        })),
    )
    .expect("leaflet operation");
}

fn add_array_operations(plan: &mut Plan, options: &PhysicalOptions, policy: &AnalysisPolicy) {
    plan.add(
        "linear",
        PlanOperation::Physical(Box::new(PhysicalRequest::LinearDensity {
            weights: 0,
            options: options.linear,
            policy: policy.clone(),
        })),
    )
    .expect("linear density operation");
    plan.add(
        "density",
        PlanOperation::Physical(Box::new(PhysicalRequest::DensityMap {
            weights: 0,
            spec: options.grid,
            policy: policy.clone(),
        })),
    )
    .expect("density operation");
    plan.add(
        "pore",
        PlanOperation::Physical(Box::new(PhysicalRequest::PoreProfile {
            radii: 0,
            options: options.pore,
            policy: policy.clone(),
        })),
    )
    .expect("pore operation");
    plan.add(
        "surface_contacts",
        PlanOperation::Physical(Box::new(PhysicalRequest::SurfaceContacts {
            radii: 0,
            tolerance: 0.5,
            probe: 1.0,
            density: 32.0,
            minimum_area: 0.2,
            backend: SpatialBackend::BruteForce,
            policy: policy.clone(),
        })),
    )
    .expect("surface contacts operation");
}

fn plan_values(
    plan: &Plan,
    structure: &crate::Structure,
    weights: &[f64],
    radii: &[f32],
) -> std::collections::BTreeMap<String, PlanValue> {
    let scalars = [ScalarInput { values: weights }];
    let floats = [FloatInput { values: radii }];
    plan.execute(
        PlanInput {
            structure: Some(structure),
            scalars: &scalars,
            floats: &floats,
            ..Default::default()
        },
        &molframe_core::ExecutionContext::default(),
    )
    .expect("physical plan")
    .entries
    .into_iter()
    .map(|entry| (entry.id.to_string(), entry.value))
    .collect()
}

fn assert_structure_matches(
    values: &std::collections::BTreeMap<String, PlanValue>,
    direct: &DirectPhysicalAnalyses,
) {
    assert!(matches!(
        values.get("radial"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::RadialDistribution(actual) if actual.value == direct.radial.value)
    ));
    assert!(matches!(
        values.get("coordination"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::CoordinationNumbers(actual) if actual.value == direct.coordination.value)
    ));
    assert!(matches!(
        values.get("leaflets"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::Leaflets(actual) if actual.value == direct.leaflets.value)
    ));
}

fn assert_array_matches(
    values: &std::collections::BTreeMap<String, PlanValue>,
    direct: &DirectPhysicalAnalyses,
) {
    assert!(matches!(
        values.get("linear"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::LinearDensity(actual) if actual.value == direct.linear.value)
    ));
    assert!(matches!(
        values.get("density"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::DensityMap(actual) if actual.value == direct.density.value)
    ));
    assert!(matches!(
        values.get("pore"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::PoreProfile(actual) if actual.value == direct.pore.value)
    ));
    assert!(matches!(
        values.get("surface_contacts"),
        Some(PlanValue::Physical(value))
            if matches!(value.as_ref(), PhysicalValue::SurfaceContacts(actual) if actual.value == direct.surface_contacts.value)
    ));
}

#[test]
fn physical_plan_matches_governed_direct_kernels() {
    let structure = structure();
    let left = crate::AtomSelection::from_sorted(vec![0, 1, 2]);
    let weights = [1.0_f64, 2.0, 1.0, 0.5];
    let radii = [1.5_f32, 1.5, 1.5, 1.5];
    let options = options();
    let policy = AnalysisPolicy::default();
    let direct = direct_analyses(&structure, &left, &weights, &radii, &options, &policy);
    let mut plan = Plan::new();
    add_structure_operations(&mut plan, &left, &options, &policy);
    add_array_operations(&mut plan, &options, &policy);
    let values = plan_values(&plan, &structure, &weights, &radii);
    assert_structure_matches(&values, &direct);
    assert_array_matches(&values, &direct);
}
