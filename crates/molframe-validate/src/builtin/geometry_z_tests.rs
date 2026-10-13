use super::*;
use crate::builtin_validation::fixtures::reference;
use molframe_bench::{Sample, structure as sample};

fn scores(which: Sample) -> Result<GeometryScores, BuiltinError> {
    let structure = sample(which);
    geometry_scores(
        &Layout::new(&structure),
        &reference(),
        4.0,
        &ExecutionContext::default(),
    )
}

#[test]
fn only_residues_with_targets_are_measured() -> Result<(), BuiltinError> {
    let s = scores(Sample::Tiny)?;
    assert_eq!(s.worst_bond_z.len(), 46);
    assert!(s.bond_count > 0 && s.angle_count > 0);
    assert!(s.bond_rmsz.is_finite() && s.angle_rmsz.is_finite());
    assert_eq!(s.tables.len(), 1);
    assert_eq!(s.tables[0].name, "synthetic-geometry");
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
    // Alanine has two bonds and one angle in-residue; every linked residue adds
    // the C-N bond and the CA-C-N angle regardless of its own name.
    assert!(s.bond_count > 60, "{}", s.bond_count);
    assert!(s.angle_count > 60, "{}", s.angle_count);
    Ok(())
}

#[test]
fn a_run_without_a_geometry_table_fails_with_the_option_to_set() {
    let structure = sample(Sample::Tiny);
    let mut set = reference();
    set.geometry = molframe_chem::TableRequest::named("not-in-the-cache");
    let result = geometry_scores(
        &Layout::new(&structure),
        &set,
        4.0,
        &ExecutionContext::default(),
    );
    assert!(matches!(result, Err(BuiltinError::Table(_))));
}

/// Parity with the published targets; runs only when the tables are present.
#[test]
#[ignore = "needs the geometry table selected by the profile or MOLFRAME_DATA_DIR"]
fn real_geometry_table_scores_crambin_tightly() -> Result<(), BuiltinError> {
    let structure = sample(Sample::Tiny);
    let s = geometry_scores(
        &Layout::new(&structure),
        &ReferenceSet::default(),
        4.0,
        &ExecutionContext::default(),
    )?;
    assert!(s.bond_rmsz < 3.0 && s.angle_rmsz < 3.0);
    Ok(())
}
