use super::*;

const GEOMETRY: &str = include_str!("../../data/geometry.tsv");
const CHI: &str = include_str!("../../data/chi_atoms.tsv");

#[test]
fn geometry_has_backbone_targets_for_every_residue() {
    let targets = GeometryTargets::parse(GEOMETRY);
    assert!(targets.len() > 300);
    let alanine = targets.residue("ALA");
    let n_ca = alanine.iter().find(|t| t.names == ["N", "CA"]);
    assert_eq!(n_ca.map(|t| t.mean), Some(1.459));
    assert!(alanine.iter().any(Target::is_angle));
}

#[test]
fn link_terms_prefer_the_proline_entries() {
    let targets = GeometryTargets::parse(GEOMETRY);
    let generic = targets.link("ALA");
    let proline = targets.link("PRO");
    assert_eq!(generic.len(), 4);
    assert_eq!(proline.len(), 4);
    let bond = |terms: &[&Target]| terms.iter().find(|t| t.names == ["C", "N"]).map(|t| t.mean);
    assert_eq!(bond(&generic), Some(1.329));
    assert_eq!(bond(&proline), Some(1.341));
}

#[test]
fn chi_definitions_cover_the_side_chains() {
    let definitions = ChiDefinitions::parse(CHI);
    let chi_definition = |name: &str| definitions.get(name);
    assert_eq!(chi_definition("ARG").map(|d| d.chi.len()), Some(4));
    assert_eq!(chi_definition("VAL").map(|d| d.chi.len()), Some(1));
    assert_eq!(
        chi_definition("TYR").map(|d| d.table.as_str()),
        Some("phetyr")
    );
    assert!(chi_definition("GLY").is_none());
}

#[test]
fn malformed_lines_are_skipped() {
    let parsed = GeometryTargets::parse("bond\tALA\tCA-CB\tx\t0.02\tEH\nshort\n");
    assert!(parsed.is_empty());
}
