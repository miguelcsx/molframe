use super::*;
use crate::builtin_validation::bundled_reference;
use molframe_bench::{Sample, structure as sample};

fn scores(which: Sample) -> Result<GeometryScores, BuiltinError> {
    let structure = sample(which);
    geometry_scores(
        &Layout::new(&structure),
        bundled_reference()?,
        4.0,
        &ExecutionContext::default(),
    )
}

#[test]
fn crambin_is_measured_with_every_term_kind() -> Result<(), BuiltinError> {
    let s = scores(Sample::Tiny)?;
    assert_eq!(s.worst_bond_z.len(), 46);
    assert!(s.bond_count > 300 && s.angle_count > 400);
    assert!(
        s.bond_rmsz.is_finite() && s.bond_rmsz < 3.0,
        "{}",
        s.bond_rmsz
    );
    assert!(
        s.angle_rmsz.is_finite() && s.angle_rmsz < 3.0,
        "{}",
        s.angle_rmsz
    );
    Ok(())
}

#[test]
fn outliers_name_their_atoms_and_exceed_the_cutoff() -> Result<(), BuiltinError> {
    let s = scores(Sample::Medium)?;
    assert!(s.outliers.iter().all(|o| o.z.abs() > 4.0));
    assert!(s.outliers.iter().all(|o| o.atoms[0] != NO_ATOM));
    assert!(s.outliers.windows(2).all(|w| w[0].residue <= w[1].residue));
    Ok(())
}

#[test]
fn link_terms_cross_the_peptide_bond() -> Result<(), BuiltinError> {
    let s = scores(Sample::Small)?;
    assert!(
        s.bond_count > 500,
        "backbone link bonds are counted: {}",
        s.bond_count
    );
    Ok(())
}
