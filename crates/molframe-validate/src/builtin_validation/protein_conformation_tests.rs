use super::*;
use crate::builtin_validation::fixtures::reference;
use molframe_bench::{Sample, structure as sample};
use molframe_core::ExecutionContext;

fn count(values: &[u8], wanted: u8) -> usize {
    values
        .iter()
        .map(|value| usize::from(*value == wanted))
        .sum()
}

fn table(which: Sample) -> Result<ProteinTable, BuiltinError> {
    let structure = sample(which);
    protein_table(
        &Layout::new(&structure),
        &reference(),
        &ExecutionContext::default(),
    )
}

#[test]
fn ubiquitin_has_no_ramachandran_outliers() -> Result<(), BuiltinError> {
    let t = table(Sample::Small)?;
    assert_eq!(count(&t.rama_class, 0), 0);
    assert!(count(&t.rama_class, 2) >= 60, "most residues are favoured");
    Ok(())
}

#[test]
fn crambin_terminal_residues_lack_one_torsion() -> Result<(), BuiltinError> {
    let t = table(Sample::Tiny)?;
    assert_eq!(t.len(), 46);
    assert!(t.phi[0].is_nan() && t.psi[0].is_finite());
    assert!(t.phi[45].is_finite() && t.psi[45].is_nan());
    assert_eq!(t.rama_class[0], NOT_EVALUATED);
    assert_eq!(t.chi.len(), 46 * 4);
    assert_eq!(count(&t.omega_class, OmegaClass::Cis.code()), 0);
    Ok(())
}

#[test]
fn only_residues_with_a_chi_definition_get_a_rotamer_class() -> Result<(), BuiltinError> {
    let structure = sample(Sample::Small);
    let layout = Layout::new(&structure);
    let t = protein_table(&layout, &reference(), &ExecutionContext::default())?;
    let valines = (0..layout.len())
        .filter(|residue| layout.name(*residue) == "VAL")
        .count();
    assert!(valines > 0);
    let evaluated = t
        .rotamer_class
        .iter()
        .filter(|c| **c != NOT_EVALUATED)
        .count();
    assert!(
        evaluated > 0 && evaluated <= valines,
        "{evaluated} of {valines}"
    );
    assert_eq!(t.tables.len(), 3);
    assert_eq!(t.tables[0].name, "synthetic-rama");
    Ok(())
}

#[test]
fn a_run_without_a_rama_table_fails_before_scoring() {
    let structure = sample(Sample::Tiny);
    let mut set = reference();
    set.rama = molframe_chem::TableRequest::named("not-in-the-cache");
    let result = protein_table(&Layout::new(&structure), &set, &ExecutionContext::default());
    assert!(matches!(result, Err(BuiltinError::Table(_))));
}

/// Parity with the wwPDB report for 1UBQ; runs only when the tables are present.
#[test]
#[ignore = "needs rama, rotamer and chi tables selected by the profile or MOLFRAME_DATA_DIR"]
fn real_tables_match_the_wwpdb_report_for_ubiquitin() -> Result<(), BuiltinError> {
    let structure = sample(Sample::Small);
    let t = protein_table(
        &Layout::new(&structure),
        &ReferenceSet::default(),
        &ExecutionContext::default(),
    )?;
    let outliers = count(&t.rotamer_class, 0);
    assert!((4..=9).contains(&outliers), "{outliers}");
    assert_eq!(count(&t.rama_class, 0), 0);
    Ok(())
}

#[test]
fn result_does_not_depend_on_worker_count() -> Result<(), Box<dyn std::error::Error>> {
    let structure = sample(Sample::Medium);
    let layout = Layout::new(&structure);
    let one = ExecutionContext::builder().worker_budget(1).build()?;
    let many = ExecutionContext::builder().worker_budget(4).build()?;
    let a = protein_table(&layout, &reference(), &one)?;
    let b = protein_table(&layout, &reference(), &many)?;
    assert_eq!(a.rama_class, b.rama_class);
    assert_eq!(a.rotamer_class, b.rotamer_class);
    assert_eq!(a.phi.len(), b.phi.len());
    Ok(())
}

#[test]
fn omega_classes_follow_the_angle() {
    assert_eq!(OmegaClass::from_degrees(179.0), OmegaClass::Trans);
    assert_eq!(OmegaClass::from_degrees(-5.0), OmegaClass::Cis);
    assert_eq!(OmegaClass::from_degrees(120.0), OmegaClass::Twisted);
}
