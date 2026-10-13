use super::*;
use molframe_core::ExecutionContext;
use molframe_core::io::{InputBuffer, ReadOptions};

const STEM: &str = include_str!("../../tests/data/1ehz_stem.cif");

fn table(policy: &NucleicPolicy) -> Result<NucleicTable, Box<dyn std::error::Error>> {
    let input = InputBuffer::from_bytes(STEM.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new())
        .map_err(|findings| format!("{findings:?}"))?;
    Ok(nucleic_table(
        &Layout::new(&structure),
        policy,
        &ExecutionContext::default(),
    )?)
}

#[test]
fn a_form_stem_is_c3_endo_with_sane_suites() -> Result<(), Box<dyn std::error::Error>> {
    let t = table(&NucleicPolicy::default())?;
    assert_eq!(t.len(), 7);
    let c3 = t
        .pucker_class
        .iter()
        .filter(|c| **c == PUCKER_C3_ENDO)
        .count();
    assert!(c3 >= 6, "{:?}", t.pucker_class);
    for phase in t.pucker_phase.iter().filter(|p| p.is_finite()) {
        assert!((0.0..360.0).contains(phase));
    }
    assert!(t.pucker_amplitude.iter().all(|a| *a > 20.0 && *a < 55.0));
    let sane = t.suite_class.iter().filter(|c| **c == 0).count();
    assert!(sane >= 5, "{:?} {:?}", t.suite_class, t.suite_flags);
    assert_eq!(
        t.suite_class[0], NOT_EVALUATED,
        "first residue has no suite"
    );
    Ok(())
}

#[test]
fn backbone_torsions_are_measured_across_links() -> Result<(), Box<dyn std::error::Error>> {
    let t = table(&NucleicPolicy::default())?;
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
    let strict = NucleicPolicy {
        c3_endo_delta: [0.0, 10.0],
        ..NucleicPolicy::default()
    };
    let t = table(&strict)?;
    let flagged = t.suite_class.iter().filter(|c| **c == 1).count();
    assert!(flagged >= 5, "{:?}", t.suite_class);
    assert!(
        t.suite_flags
            .iter()
            .any(|f| f & SUITE_DELTA != 0 && *f != NOT_EVALUATED)
    );
    assert_eq!(t.tables, vec![strict.info]);
    Ok(())
}

#[test]
fn proteins_are_not_nucleotides() -> Result<(), Box<dyn std::error::Error>> {
    let structure = molframe_bench::structure(molframe_bench::Sample::Tiny);
    let t = nucleic_table(
        &Layout::new(&structure),
        &NucleicPolicy::default(),
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
