use crate::{
    AnalysisPolicy, BondOrder, Component, ComponentAtom, ComponentBond, ComponentKind,
    DictionaryVersion, Element, MemoryProvider, Plan, PlanInput, PlanOperation, PlanValue,
    PolymerLinkPolicy, ReadOptions, StructureRequest,
};
use std::sync::Arc;

const GNM_CIF: &str = "data_gnm\n\
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
ATOM 1 C CA . GLY A 1 1 ? 0 0 0 1.0 10.0 1 GLY A CA 1\n\
ATOM 2 C CB . GLY A 1 1 ? 1 0 0 1.0 10.0 1 GLY A CB 1\n\
ATOM 3 C CG . GLY A 1 1 ? 2 0 0 1.0 10.0 1 GLY A CG 1\n";

#[test]
fn gaussian_network_plan_matches_direct_kernel() {
    let (structure, _) = crate::read_bytes(
        GNM_CIF.as_bytes().to_vec(),
        Some("gnm.cif"),
        &ReadOptions::new(),
    )
    .expect("GNM fixture");
    let sites = crate::AtomSelection::from_sorted(vec![0, 1, 2]);
    let options = molframe_analysis::GnmOptions {
        contact_distance: 1.1,
        mode_count: 2,
        zero_mode_tolerance: 1e-10,
        memory_limit_bytes: 1_024,
        backend: crate::SpatialBackend::BruteForce,
        reduction: molframe_core::parallel::ReductionPolicy::Deterministic,
    };
    let direct = molframe_analysis::gaussian_network_model(
        structure.positions(),
        &sites,
        options,
        None,
        &molframe_core::ExecutionContext::default(),
    )
    .expect("direct GNM");
    let mut plan = Plan::new();
    plan.add(
        "gnm",
        PlanOperation::Structure(Box::new(StructureRequest::GaussianNetworkModel {
            sites,
            options,
            periodic: false,
            policy: AnalysisPolicy::default(),
        })),
    )
    .expect("GNM operation");
    let result = plan
        .execute(
            PlanInput {
                structure: Some(&structure),
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("GNM plan");
    let Some(PlanValue::Structure(value)) = result.entries.first().map(|entry| &entry.value) else {
        panic!("missing GNM result");
    };
    let super::super::StructureValue::GaussianNetworkModel(analysis) = value.as_ref() else {
        panic!("unexpected GNM result type");
    };
    assert_eq!(analysis.value, direct);
    assert_eq!(analysis.status, crate::Status::Complete);
}

const BASE_PAIR_CIF: &str = "data_base_pair\n\
loop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
ATOM 1 N D ADE A 1 0 0 0\n\
ATOM 2 H H ADE A 1 1 0 0\n\
ATOM 3 O A URA B 1 2.8 0 0\n";

#[test]
fn base_pair_plan_matches_direct_kernel_with_an_explicit_ccd_provider() {
    let (parsed, _) = crate::read_bytes(
        BASE_PAIR_CIF.as_bytes().to_vec(),
        Some("base-pair.cif"),
        &ReadOptions::new(),
    )
    .expect("base-pair fixture");
    let provider: Arc<dyn crate::ComponentProvider> = Arc::new(
        MemoryProvider::new(
            DictionaryVersion::new("test"),
            [base_pair_donor_component(), base_pair_acceptor_component()],
        )
        .expect("component fixtures are unique"),
    );
    let structure =
        crate::apply_component_chemistry(&parsed, provider.as_ref(), PolymerLinkPolicy::Disabled)
            .expect("component chemistry")
            .structure;
    let options = crate::analysis::BasePairOptions {
        hydrogen_bonds: crate::analysis::HydrogenBondOptions {
            maximum_donor_acceptor_distance: 3.5,
            minimum_angle_degrees: 150.0,
            backend: crate::SpatialBackend::BruteForce,
            periodic: false,
        },
        minimum_hydrogen_bonds: 1,
    };
    let direct = crate::analysis::base_pairs(
        &structure,
        provider.as_ref(),
        options,
        &molframe_core::ExecutionContext::default(),
    )
    .expect("direct base-pair kernel");
    let mut plan = Plan::new();
    plan.add(
        "base_pairs",
        PlanOperation::Structure(Box::new(StructureRequest::BasePairs {
            provider,
            options,
            policy: AnalysisPolicy::default(),
        })),
    )
    .expect("base-pair operation");
    let result = plan
        .execute(
            PlanInput {
                structure: Some(&structure),
                ..Default::default()
            },
            &molframe_core::ExecutionContext::default(),
        )
        .expect("base-pair plan");
    let Some(PlanValue::Structure(value)) = result.entries.first().map(|entry| &entry.value) else {
        panic!("missing base-pair result");
    };
    let super::super::StructureValue::BasePairs(analysis) = value.as_ref() else {
        panic!("unexpected base-pair result type");
    };
    assert_eq!(analysis.value, direct);
    assert_eq!(analysis.value.len(), 1);
    assert_eq!(analysis.status, crate::Status::Complete);
}

fn base_pair_donor_component() -> Component {
    Component {
        id: "ADE".into(),
        name: "adenine test".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(b'A'),
        formula: None,
        atoms: Arc::from([
            base_pair_atom("D", Element::NITROGEN),
            base_pair_atom("H", Element::HYDROGEN),
        ]),
        bonds: Arc::from([ComponentBond {
            atom_a: "D".into(),
            atom_b: "H".into(),
            order: BondOrder::Single,
            aromatic: false,
            stereo: None,
        }]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn base_pair_acceptor_component() -> Component {
    Component {
        id: "URA".into(),
        name: "uracil test".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(b'U'),
        formula: None,
        atoms: Arc::from([base_pair_atom("A", Element::OXYGEN)]),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn base_pair_atom(name: &str, element: Element) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: None,
    }
}
