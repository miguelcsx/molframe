use super::*;
use crate::builtin_validation::BuiltinError;
use crate::builtin_validation::bundled_reference;

#[test]
fn rama_helix_is_favoured_and_forbidden_corner_is_not() -> Result<(), BuiltinError> {
    let grids = &bundled_reference()?.rama;
    assert!(grids.density(RamaCategory::General, -63.0, -43.0) >= RAMA_FAVOURED);
    assert!(grids.density(RamaCategory::General, -63.0, -43.0) > 0.0);
    assert!(
        grids.density(RamaCategory::General, -179.0, 1.0) < RamaCategory::General.outlier_density()
    );
    assert!(grids.density(RamaCategory::Gly, 80.0, 0.0) >= RAMA_FAVOURED);
    Ok(())
}

#[test]
fn rama_wraps_angles() -> Result<(), BuiltinError> {
    let grids = &bundled_reference()?.rama;
    let wrapped = grids.density(RamaCategory::General, 297.0, 317.0);
    assert!((grids.density(RamaCategory::General, -63.0, -43.0) - wrapped).abs() < 1e-6);
    Ok(())
}

#[test]
fn rotamer_tables_have_all_seventeen() -> Result<(), BuiltinError> {
    let tables = &bundled_reference()?.rotamers;
    assert_eq!(tables.len(), 17);
    assert_eq!(tables.get("arg").map(RotamerTable::dimensions), Some(4));
    assert_eq!(tables.get("phetyr").map(RotamerTable::dimensions), Some(2));
    assert_eq!(tables.get("val").map(RotamerTable::dimensions), Some(1));
    Ok(())
}

#[test]
fn valine_trans_is_favoured_and_eclipsed_is_not() -> Result<(), BuiltinError> {
    let valine = bundled_reference()?.rotamers.get("val");
    assert_eq!(
        valine.and_then(|t| t.classify(&[175.0])),
        Some(RotamerClass::Favoured)
    );
    assert_eq!(
        valine.and_then(|t| t.classify(&[0.0])),
        Some(RotamerClass::Outlier)
    );
    assert_eq!(valine.and_then(|t| t.classify(&[])), None);
    Ok(())
}

#[test]
fn malformed_tables_are_rejected() {
    assert!(RamaGrids::from_bytes(b"nope").is_err());
    assert!(RamaGrids::from_bytes(b"MFRM\x06").is_err());
    assert!(RotamerTables::from_bytes(b"MFRT\x01").is_err());
}

#[test]
fn thresholds_decide_the_class_and_can_be_replaced() -> Result<(), BuiltinError> {
    let default = bundled_reference()?.rama.thresholds();
    assert!((default.favoured - RAMA_FAVOURED).abs() < f32::EPSILON);
    let grids = RamaGrids::from_bytes(include_bytes!("../../data/rama.bin"))?;
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
