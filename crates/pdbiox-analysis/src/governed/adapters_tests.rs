use super::{
    GovernedAnalysisError, GovernedNativeError, PhysicalKernelError, analyse_structure,
    contacts_kernel, linear_density_kernel, native_contact_fraction_kernel,
};
use crate::{CartesianAxis, LinearDensityOptions};
use pdbiox_core::ExecutionContext;
use pdbiox_core::contract::{AnalysisPolicy, PeriodicPolicy};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_spatial::SpatialBackend;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1 0 0\n";

#[test]
fn concrete_contact_adapter_records_identity_and_parameters() {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let Ok((structure, _)) = pdbiox_cif::read(&input, &ReadOptions::new()) else {
        panic!("valid structure fixture");
    };
    let kernel = contacts_kernel(1.5, SpatialBackend::BruteForce);
    let Ok(result) = analyse_structure(
        &structure,
        &AnalysisPolicy::default(),
        &kernel,
        &ExecutionContext::default(),
    ) else {
        panic!("contact analysis should succeed");
    };
    assert_eq!(result.value.len(), 1);
    assert_eq!(result.coverage.used, 2);
    assert_eq!(
        result
            .provenance
            .algorithm
            .as_ref()
            .map(pdbiox_core::contract::AlgorithmId::name),
        Some("atom-contacts")
    );
    assert!(result.provenance.parameters.contains_key("cutoff"));
    assert!(result.provenance.parameters.contains_key("spatial_backend"));
}

fn structure() -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let Ok((structure, _)) = pdbiox_cif::read(&input, &ReadOptions::new()) else {
        panic!("valid structure fixture");
    };
    structure
}

#[test]
fn atom_aligned_side_inputs_are_never_padded() {
    let structure = structure();
    let kernel = linear_density_kernel(
        &[1.0],
        LinearDensityOptions {
            axis: CartesianAxis::X,
            minimum: 0.0,
            maximum: 2.0,
            bins: 2,
        },
    );
    let error = analyse_structure(
        &structure,
        &AnalysisPolicy::default(),
        &kernel,
        &ExecutionContext::default(),
    );
    assert!(matches!(
        error,
        Err(GovernedAnalysisError::Kernel(
            PhysicalKernelError::SideInputLength
        ))
    ));
}

#[test]
fn native_q_rejects_unimplemented_periodic_semantics() {
    let structure = structure();
    let kernel = native_contact_fraction_kernel(&structure, 2.0, 1.0, SpatialBackend::BruteForce);
    let policy = AnalysisPolicy {
        periodic: PeriodicPolicy::MinimumImage,
        ..AnalysisPolicy::default()
    };
    let error = analyse_structure(&structure, &policy, &kernel, &ExecutionContext::default());
    assert!(matches!(
        error,
        Err(GovernedAnalysisError::Kernel(
            GovernedNativeError::PeriodicUnsupported
        ))
    ));
}
