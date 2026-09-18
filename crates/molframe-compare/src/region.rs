//! Typed RMSD results for explicitly mapped structural regions.

use crate::{ComparisonAlignment, DistanceMeasurement, PointMapping, measure_mapping};

/// RMSD over an explicitly mapped interface region.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InterfaceRmsd(pub f64);

/// RMSD over an explicitly mapped binding-pocket region.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PocketRmsd(pub f64);

fn region_measurement(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    alignment: ComparisonAlignment,
) -> DistanceMeasurement {
    measure_mapping(reference, model, mapping, alignment)
}

/// Measures an explicitly mapped interface after the caller-selected alignment.
#[must_use]
pub fn interface_rmsd(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    alignment: ComparisonAlignment,
) -> InterfaceRmsd {
    InterfaceRmsd(region_measurement(reference, model, mapping, alignment).rmsd)
}

/// Measures an explicitly mapped binding pocket after caller-selected alignment.
#[must_use]
pub fn pocket_rmsd(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    alignment: ComparisonAlignment,
) -> PocketRmsd {
    PocketRmsd(region_measurement(reference, model, mapping, alignment).rmsd)
}

#[cfg(test)]
#[path = "region_tests.rs"]
mod tests;
