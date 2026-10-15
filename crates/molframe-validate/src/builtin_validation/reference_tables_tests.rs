use super::*;
use crate::builtin_validation::fixtures::{CHI, GEOMETRY, identity};

fn geometry() -> Result<Provided<GeometryTargets>, TableTextError> {
    GeometryTargets::decode(GEOMETRY.as_bytes(), identity("fallback"))
}

fn reject(text: &str) -> bool {
    GeometryTargets::decode(text.as_bytes(), identity("fallback")).is_err()
}

#[test]
fn geometry_rows_land_under_their_residue() -> Result<(), TableTextError> {
    let provided = geometry()?;
    assert_eq!(provided.identity.name, "synthetic-geometry");
    let alanine = provided.table.residue("ALA");
    let n_ca = alanine.iter().find(|t| t.names == ["N", "CA"]);
    assert_eq!(n_ca.map(|t| t.mean), Some(1.46));
    assert!(alanine.iter().any(Target::is_angle));
    assert!(provided.table.residue("GLY").is_empty());
    assert_eq!(provided.table.len(), 6);
    Ok(())
}

#[test]
fn link_terms_prefer_the_entries_of_the_following_residue() -> Result<(), TableTextError> {
    let targets = geometry()?.table;
    let generic: Vec<&Target> = targets.link("ALA").collect();
    let special: Vec<&Target> = targets.link("PRO").collect();
    assert_eq!(generic.len(), 2);
    assert_eq!(special.len(), 2);
    let bond = |terms: &[&Target]| terms.iter().find(|t| t.names == ["C", "N"]).map(|t| t.mean);
    assert_eq!(bond(&generic), Some(1.33));
    assert_eq!(bond(&special), Some(1.34));
    Ok(())
}

#[test]
fn chi_definitions_name_their_rotamer_table() -> Result<(), TableTextError> {
    let definitions = ChiDefinitions::decode(CHI.as_bytes(), identity("fallback"))?.table;
    assert_eq!(definitions.get("VAL").map(|d| d.chi.len()), Some(1));
    assert_eq!(
        definitions.get("VAL").map(|d| d.table.as_str()),
        Some("val")
    );
    assert!(definitions.get("GLY").is_none());
    assert_eq!(definitions.len(), 1);
    Ok(())
}

#[test]
fn malformed_rows_are_errors_not_skipped() {
    let head = "name: x\nversion: 1\nsource: s\nlicence: CC0-1.0\n";
    assert!(reject(&format!("{head}bond\tALA\tCA-CB\tx\t0.02\n")));
    assert!(reject(&format!("{head}short\n")));
    assert!(reject(&format!("{head}bond\tALA\tCA-CB\t1.5\t0\n")));
    assert!(reject(&format!("{head}bond\tALA\tCA-CB-CG\t1.5\t0.02\n")));
    assert!(reject(&format!("{head}bond\t*\tCA-CB\t1.5\t0.02\n")));
    assert!(reject(head));
    let bad_chi = format!("{head}VAL\tval\tN CA CB\n");
    assert!(ChiDefinitions::decode(bad_chi.as_bytes(), identity("f")).is_err());
    let twice = format!("{head}VAL\tval\tN CA CB CG1\nVAL\tval\tN CA CB CG2\n");
    assert!(ChiDefinitions::decode(twice.as_bytes(), identity("f")).is_err());
}
