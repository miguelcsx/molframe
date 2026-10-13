use super::*;
use crate::builtin_validation::fixtures::{identity, rama_bytes, rotamer_bytes};

fn rama() -> Result<Provided<RamaGrids>, TableTextError> {
    RamaGrids::decode(&rama_bytes(), identity("fallback"))
}

fn rotamers() -> Result<Provided<RotamerTables>, TableTextError> {
    RotamerTables::decode(&rotamer_bytes(), identity("fallback"))
}

#[test]
fn rama_density_follows_the_grid_and_the_forbidden_cell_is_an_outlier() -> Result<(), TableTextError>
{
    let grids = rama()?.table;
    assert!(grids.density(RamaCategory::General, -63.0, -43.0) > 0.99);
    let forbidden = grids.density(RamaCategory::General, -179.0, 1.0);
    assert!(forbidden.abs() < f32::EPSILON);
    assert_eq!(grids.classify(RamaCategory::General, forbidden), 0);
    assert_eq!(grids.classify(RamaCategory::Gly, 1.0), 2);
    Ok(())
}

#[test]
fn rama_wraps_angles() -> Result<(), TableTextError> {
    let grids = rama()?.table;
    let wrapped = grids.density(RamaCategory::General, 297.0, 317.0);
    assert!((grids.density(RamaCategory::General, -63.0, -43.0) - wrapped).abs() < 1e-6);
    Ok(())
}

#[test]
fn rama_identity_and_thresholds_come_from_the_header() -> Result<(), TableTextError> {
    let provided = rama()?;
    assert_eq!(provided.identity, identity("synthetic-rama"));
    let cuts = provided.table.thresholds();
    assert!((cuts.favoured - 0.02).abs() < f32::EPSILON);
    assert!((cuts.outlier[3] - 0.002).abs() < f32::EPSILON);
    Ok(())
}

#[test]
fn thresholds_decide_the_class_and_can_be_replaced() -> Result<(), TableTextError> {
    let grids = rama()?.table;
    let default = grids.thresholds();
    assert_eq!(grids.classify(RamaCategory::General, 0.0001), 0);
    assert_eq!(grids.classify(RamaCategory::General, 0.001), 1);
    assert_eq!(grids.classify(RamaCategory::General, 0.05), 2);
    let strict = grids.with_thresholds(RamaThresholds {
        favoured: 0.5,
        ..default
    });
    assert_eq!(strict.classify(RamaCategory::General, 0.05), 1);
    Ok(())
}

#[test]
fn a_rama_table_without_cut_offs_is_refused() {
    let cells = RAMA_BINS * RAMA_BINS;
    let mut payload = vec![6u8];
    for _ in 0..6 {
        payload.extend_from_slice(&[0u8; 8]);
        payload.extend(std::iter::repeat_n(255u8, cells));
    }
    let bare =
        molframe_chem::join_binary(*b"MFR2", &identity("bare"), "", &payload).unwrap_or_default();
    assert!(RamaGrids::decode(&bare, identity("fallback")).is_err());
}

#[test]
fn rotamer_tables_are_found_by_name() -> Result<(), TableTextError> {
    let tables = rotamers()?.table;
    assert_eq!(tables.len(), 1);
    assert_eq!(tables.get("val").map(RotamerTable::dimensions), Some(1));
    assert!(tables.get("arg").is_none());
    Ok(())
}

#[test]
fn rotamer_classes_follow_the_cells() -> Result<(), TableTextError> {
    let tables = rotamers()?.table;
    let valine = tables.get("val");
    assert_eq!(
        valine.and_then(|t| t.classify(&[175.0])),
        Some(RotamerClass::Favoured)
    );
    assert_eq!(
        valine.and_then(|t| t.classify(&[0.0])),
        Some(RotamerClass::Outlier)
    );
    assert_eq!(
        valine.and_then(|t| t.classify(&[300.0])),
        Some(RotamerClass::Allowed)
    );
    assert_eq!(valine.and_then(|t| t.classify(&[])), None);
    Ok(())
}

#[test]
fn malformed_tables_are_rejected() {
    let fallback = || identity("fallback");
    assert!(RamaGrids::decode(b"nope", fallback()).is_err());
    assert!(RamaGrids::decode(b"MFR2\x06", fallback()).is_err());
    assert!(RotamerTables::decode(b"MFT2\x01", fallback()).is_err());
    let wrong_magic = rotamer_bytes();
    assert!(RamaGrids::decode(&wrong_magic, fallback()).is_err());
}
