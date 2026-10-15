use super::sasa_kernel;
use crate::{GovernedAnalysisError, analyse_structure};
use molframe_chem::{RadiusSelection, TableRequest, TableStore, VdwRadii};
use molframe_cif::read;
use molframe_core::contract::{Analysis, AnalysisPolicy, PolicyField, RadiiSet};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{ExecutionContext, InputBuffer, ReadOptions, Structure};
use molframe_surface::SasaOptions;

const TEXT: &str = "data_t\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 GLY A 1 0 0 0\nATOM 2 C C2 GLY A 2 3.0 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(TEXT.as_bytes().to_vec());
    match read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn policy(radii: RadiiSet) -> AnalysisPolicy {
    AnalysisPolicy {
        vdw_radii: radii,
        ..AnalysisPolicy::default()
    }
}

fn options(request: &TableRequest<VdwRadii>) -> SasaOptions {
    match SasaOptions::from_words(
        None,
        None,
        None,
        Some("all"),
        None,
        request,
        molframe_bench::table_store(),
    ) {
        Ok(options) => options,
        Err(error) => panic!("the options are valid: {error}"),
    }
}

/// Runs the kernel the way a front end does: the selection chosen from the arguments,
/// resolved against the policy, then handed to the kernel.
fn run(
    argument: Option<&str>,
    policy_supplied: bool,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<f64>>, String> {
    let selection = RadiusSelection::choose(argument, policy_supplied);
    let store: &TableStore = molframe_bench::table_store();
    let request = selection
        .request(policy.vdw_radii, store)
        .map_err(|error| Diagnostic::from(error).code().to_string())?;
    let kernel = sasa_kernel(options(&request), selection);
    match analyse_structure(&structure(), policy, &kernel, &ExecutionContext::default()) {
        Ok(result) => Ok(result),
        Err(GovernedAnalysisError::Kernel(error)) => Err(error.to_string()),
        Err(error) => panic!("unexpected failure: {error}"),
    }
}

fn areas(result: &Analysis<Vec<f64>>) -> Vec<f64> {
    match result.value() {
        Some(areas) => areas.clone(),
        None => panic!("the analysis has a value"),
    }
}

fn radii_assumption(result: &Analysis<Vec<f64>>) -> String {
    match result
        .assumptions
        .iter()
        .find(|assumption| assumption.field == PolicyField::VdwRadii)
    {
        Some(assumption) => assumption.value.to_string(),
        None => panic!("the result carries a vdw_radii assumption"),
    }
}

#[test]
fn the_policy_alone_picks_the_table_and_the_areas_differ_by_set() {
    let bondi = run(None, true, &policy(RadiiSet::Bondi));
    let charmm = run(None, true, &policy(RadiiSet::Charmm));
    let (Ok(bondi), Ok(charmm)) = (bondi, charmm) else {
        panic!("both sets are in the data directory");
    };
    assert_ne!(areas(&bondi), areas(&charmm));
    assert!(radii_assumption(&bondi).contains("bondi"));
    assert!(radii_assumption(&charmm).contains("charmm"));
    assert!(!radii_assumption(&charmm).contains("bondi"));
}

#[test]
fn the_policy_table_equals_the_same_table_asked_for_directly() {
    let by_policy = run(None, true, &policy(RadiiSet::Charmm));
    let by_argument = run(Some("charmm"), false, &AnalysisPolicy::default());
    let (Ok(by_policy), Ok(by_argument)) = (by_policy, by_argument) else {
        panic!("both resolve");
    };
    assert_eq!(areas(&by_policy), areas(&by_argument));
}

#[test]
fn an_argument_equal_to_the_policy_is_accepted() {
    let agreed = run(Some("charmm"), true, &policy(RadiiSet::Charmm));
    let Ok(agreed) = agreed else {
        panic!("an argument naming the policy's table is not a conflict");
    };
    assert!(radii_assumption(&agreed).contains("charmm"));
}

#[test]
fn an_argument_that_names_another_table_than_the_policy_is_a_typed_error() {
    let refused = run(Some("bondi"), true, &policy(RadiiSet::Charmm));
    assert_eq!(refused.err(), Some(Code::E7103.to_string()));
}

#[test]
fn a_table_resolved_under_another_policy_is_refused_when_the_kernel_runs() {
    let request = TableRequest::named("bondi");
    let kernel = sasa_kernel(options(&request), RadiusSelection::Policy);
    let result = analyse_structure(
        &structure(),
        &policy(RadiiSet::Charmm),
        &kernel,
        &ExecutionContext::default(),
    );
    assert!(matches!(result, Err(GovernedAnalysisError::Kernel(_))));
}

#[test]
fn a_policy_less_call_keeps_the_profile_table() {
    let result = run(None, false, &AnalysisPolicy::default());
    let Ok(result) = result else {
        panic!("the profile supplies the table");
    };
    assert!(!radii_assumption(&result).is_empty());
}
