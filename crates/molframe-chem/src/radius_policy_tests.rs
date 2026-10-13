use super::*;
use molframe_bench::{Sample, structure};
use molframe_core::Element;

fn finite(radii: &[f32]) -> usize {
    radii.iter().filter(|radius| radius.is_finite()).count()
}

fn count_where(structure: &Structure, keep: impl Fn(&str) -> bool) -> usize {
    structure
        .data()
        .residues()
        .filter(|residue| residue.name().is_some_and(&keep))
        .map(|residue| residue.atoms().count())
        .sum()
}

#[test]
fn default_measures_the_polymer_only() {
    let hemoglobin = structure(Sample::Medium);
    let radii = policy_radii(&hemoglobin, &RadiusPolicy::default());
    let hetero = count_where(&hemoglobin, |name| name == "HEM" || name == "HOH")
        + count_where(&hemoglobin, |name| name == "PO4");
    assert_eq!(finite(&radii), hemoglobin.atom_count() as usize - hetero);
    assert!(hetero > 0);
}

#[test]
fn classes_switch_independently() {
    let hemoglobin = structure(Sample::Medium);
    let total = hemoglobin.atom_count() as usize;
    let waters = count_where(&hemoglobin, is_water_component);
    assert_eq!(waters, 221);
    let policy = |hetero| RadiusPolicy {
        hetero,
        metals: MetalFallback::Vdw,
        ..RadiusPolicy::default()
    };
    let all = policy_radii(&hemoglobin, &policy(HeteroPolicy::ALL));
    assert_eq!(finite(&all), total);
    let solute = policy_radii(&hemoglobin, &policy(HeteroPolicy::SOLUTE));
    assert_eq!(finite(&solute), total - waters);
}

#[test]
fn metals_follow_the_stated_fallback() {
    let hemoglobin = structure(Sample::Medium);
    let iron = hemoglobin
        .data()
        .atoms()
        .find(|atom| atom.element() == Some(Element::IRON))
        .map(|atom| atom.index().as_usize())
        .expect("4HHB has haem iron");
    let radius = |metals| {
        let policy = RadiusPolicy {
            metals,
            hetero: HeteroPolicy::SOLUTE,
            ..RadiusPolicy::default()
        };
        policy_radii(&hemoglobin, &policy)[iron]
    };
    assert!((radius(MetalFallback::Vdw) - 2.44).abs() < 1e-6);
    assert!((radius(MetalFallback::Ionic) - 0.78).abs() < 1e-4);
    assert!(radius(MetalFallback::None).is_nan());
}

#[test]
fn a_set_value_always_beats_the_fallback() {
    let crambin = structure(Sample::Tiny);
    let radii = policy_radii(&crambin, &RadiusPolicy::default());
    let bondi = crate::atom_radii(&crambin, RadiusSet::Bondi);
    for (policy, plain) in radii.iter().zip(&bondi) {
        assert_eq!(policy.to_bits(), plain.to_bits());
    }
}

#[test]
fn modified_residues_written_as_heterogens_stay_polymer() {
    assert!(HeteroPolicy::POLYMER.includes(ResidueClass::Polymer));
    assert!(!HeteroPolicy::POLYMER.includes(ResidueClass::Water));
    assert!(HeteroPolicy::ALL.includes(ResidueClass::Ion));
}
