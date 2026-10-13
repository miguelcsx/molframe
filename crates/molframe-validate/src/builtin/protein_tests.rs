use super::*;
use crate::builtin_validation::bundled_reference;
use molframe_bench::{Sample, structure as sample};
use molframe_core::ExecutionContext;

fn count(values: &[u8], wanted: u8) -> usize {
    values.iter().filter(|value| **value == wanted).count()
}

fn table(which: Sample) -> Result<ProteinTable, BuiltinError> {
    let structure = sample(which);
    protein_table(
        &Layout::new(&structure),
        bundled_reference()?,
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
    assert!(count(&t.omega_class, OmegaClass::Cis.code()) == 0);
    Ok(())
}

#[test]
fn ubiquitin_rotamer_outliers_match_the_wwpdb_report() -> Result<(), BuiltinError> {
    let t = table(Sample::Small)?;
    let evaluated = t
        .rotamer_class
        .iter()
        .filter(|c| **c != NOT_EVALUATED)
        .count();
    let outliers = count(&t.rotamer_class, 0);
    assert!(evaluated >= 40);
    // The wwPDB report for 1UBQ lists 6 rotamer outliers among 68 residues.
    assert!((4..=9).contains(&outliers), "{outliers} of {evaluated}");
    Ok(())
}

#[test]
fn result_does_not_depend_on_worker_count() -> Result<(), Box<dyn std::error::Error>> {
    let structure = sample(Sample::Medium);
    let layout = Layout::new(&structure);
    let one = ExecutionContext::builder().worker_budget(1).build()?;
    let many = ExecutionContext::builder().worker_budget(4).build()?;
    let a = protein_table(&layout, bundled_reference()?, &one)?;
    let b = protein_table(&layout, bundled_reference()?, &many)?;
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
