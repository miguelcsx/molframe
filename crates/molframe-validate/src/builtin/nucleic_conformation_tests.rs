use super::*;
use crate::builtin_validation::fixtures::{nucleic_with_delta, reference};
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};

const STEM: &str = include_str!("../../tests/data/1ehz_stem.cif");

fn table(reference: &ReferenceSet) -> Result<NucleicTable, Box<dyn std::error::Error>> {
    let input = InputBuffer::from_bytes(STEM.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new())
        .map_err(|findings| format!("{findings:?}"))?;
    Ok(nucleic_table(
        &Layout::new(&structure),
        reference,
        &ExecutionContext::default(),
    )?)
}

#[test]
fn a_form_stem_is_c3_endo_with_sane_suites() -> Result<(), Box<dyn std::error::Error>> {
    let t = table(&reference())?;
    assert_eq!(t.len(), 7);
    let c3 = t
        .pucker_class
        .iter()
        .map(|c| usize::from(*c == PUCKER_C3_ENDO))
        .sum::<usize>();
    assert!(c3 >= 6, "{:?}", t.pucker_class);
    for phase in t.pucker_phase.iter().filter(|p| p.is_finite()) {
        assert!((0.0..360.0).contains(phase));
    }
    assert!(t.pucker_amplitude.iter().all(|a| *a > 20.0 && *a < 55.0));
    let sane = t
        .suite_class
        .iter()
        .map(|c| usize::from(*c == 0))
        .sum::<usize>();
    assert!(sane >= 5, "{:?} {:?}", t.suite_class, t.suite_flags);
    assert_eq!(
        t.suite_class[0], NOT_EVALUATED,
        "first residue has no suite"
    );
    Ok(())
}

#[test]
fn backbone_torsions_are_measured_across_links() -> Result<(), Box<dyn std::error::Error>> {
    let t = table(&reference())?;
    let at = |residue: usize, which: usize| t.torsions[residue * NUCLEIC_TORSIONS + which];
    assert!(at(0, 0).is_nan(), "no alpha without a preceding phosphate");
    assert!(at(3, 0).is_finite() && at(3, 4).is_finite() && at(3, 5).is_finite());
    assert!(at(6, 4).is_nan(), "no epsilon after the last residue");
    // A-form: gamma gauche+, delta about 80, chi anti near -160.
    assert!((30.0..90.0).contains(&at(3, 2)), "{}", at(3, 2));
    assert!((55.0..110.0).contains(&at(3, 3)), "{}", at(3, 3));
    assert!(at(3, 6) < -140.0 || at(3, 6) > 170.0, "{}", at(3, 6));
    Ok(())
}

#[test]
fn a_stricter_policy_changes_the_verdict_and_is_reported() -> Result<(), Box<dyn std::error::Error>>
{
    let mut strict = reference();
    strict.nucleic = nucleic_with_delta([0.0, 10.0]);
    let t = table(&strict)?;
    let flagged = t
        .suite_class
        .iter()
        .map(|c| usize::from(*c == 1))
        .sum::<usize>();
    assert!(flagged >= 5, "{:?}", t.suite_class);
    assert!(
        t.suite_flags
            .iter()
            .any(|f| f & SUITE_DELTA != 0 && *f != NOT_EVALUATED)
    );
    assert_eq!(t.tables.len(), 1);
    assert_eq!(t.tables[0].name, "strict-nucleic");
    Ok(())
}

#[test]
fn proteins_are_not_nucleotides() -> Result<(), Box<dyn std::error::Error>> {
    let structure = molframe_bench::structure(molframe_bench::Sample::Tiny);
    let t = nucleic_table(
        &Layout::new(&structure),
        &reference(),
        &ExecutionContext::default(),
    )?;
    assert!(t.pucker_class.iter().all(|c| *c == NOT_EVALUATED));
    assert!(t.pucker_phase.iter().all(|p| p.is_nan()));
    Ok(())
}

#[test]
fn windows_wrap_around_zero() {
    assert!(within(350.0, [-18.0, 54.0]));
    assert!(within(10.0, [-18.0, 54.0]));
    assert!(!within(100.0, [-18.0, 54.0]));
}

#[test]
fn a_run_without_a_policy_table_fails_with_the_option_to_set() {
    let mut set = reference();
    set.nucleic = molframe_chem::TableRequest::named("not-in-the-cache");
    let result = table(&set);
    assert!(result.is_err());
}
