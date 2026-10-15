use super::{ResidueKernelError, residue_sasa_kernel};
use crate::cif_fixture::read_cif;
use crate::{GovernedAnalysisError, analyse_structure};
use molframe_chem::{RadiusSelection, TableRequest, VdwRadii};
use molframe_core::contract::{AnalysisPolicy, PolicyField, RadiiSet};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{ExecutionContext, Structure};
use molframe_surface::SasaOptions;

const TEXT: &str = "data_s\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N ALA A 1 0 0 0\nATOM 2 C CA ALA A 1 1.5 0 0\n\
ATOM 3 C C GLY A 2 3.0 0 0\n";

fn structure() -> Structure {
    read_cif(TEXT)
}

fn run(
    table: Option<&str>,
    supplied: bool,
    policy: &AnalysisPolicy,
) -> Result<molframe_core::contract::Analysis<crate::ResidueAreas>, Code> {
    let selection = RadiusSelection::choose(table, supplied);
    let request: TableRequest<VdwRadii> =
        match selection.request(policy.vdw_radii, molframe_bench::table_store()) {
            Ok(request) => request,
            Err(error) => return Err(Diagnostic::from(error).code()),
        };
    let options = match SasaOptions::resolve(&request, molframe_bench::table_store()) {
        Ok(options) => options,
        Err(error) => return Err(Diagnostic::from(error).code()),
    };
    let kernel = residue_sasa_kernel(options, selection);
    match analyse_structure(&structure(), policy, &kernel, &ExecutionContext::default()) {
        Ok(result) => Ok(result),
        Err(GovernedAnalysisError::Kernel(error)) => Err(Diagnostic::from(error).code()),
        Err(error) => panic!("unexpected failure: {error}"),
    }
}

fn policy(set: RadiiSet) -> AnalysisPolicy {
    AnalysisPolicy {
        vdw_radii: set,
        ..AnalysisPolicy::default()
    }
}

#[test]
fn residues_are_summed_and_the_applied_table_is_recorded() {
    let Ok(result) = run(None, true, &policy(RadiiSet::Bondi)) else {
        panic!("the bondi set is in the data directory");
    };
    let Some(areas) = result.value() else {
        panic!("a determinate result");
    };
    assert_eq!(areas.residues.len(), 2);
    assert_eq!(areas.residues[0].atoms, 2);
    let applied = result
        .assumptions
        .iter()
        .find(|assumption| assumption.field == PolicyField::VdwRadii);
    assert!(applied.is_some_and(|a| a.value.to_string().contains("bondi")));
}

#[test]
fn a_table_that_conflicts_with_the_policy_is_a_typed_error() {
    let refused = run(Some("bondi"), true, &policy(RadiiSet::Charmm));
    assert_eq!(refused.err(), Some(Code::E7103));
}

#[test]
fn the_error_type_classifies_a_missing_atom_set() {
    let error = ResidueKernelError::Residues(crate::ResidueSasaError::NoAtoms);
    assert_eq!(Diagnostic::from(error).code(), Code::E5103);
}
