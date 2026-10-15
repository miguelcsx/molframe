//! AMBER `NetCDF` schema vocabulary and unit validation.

use netcdf_reader::{NcAttrValue, NcFile, NcVariable};

use super::AmberNetcdfError;

pub(crate) const CONVENTION_ATTRIBUTE: &str = "Conventions";
pub(crate) const CONVENTION: &str = "AMBER";
pub(crate) const VERSION_ATTRIBUTE: &str = "ConventionVersion";
pub(crate) const VERSION: &str = "1.0";
pub(crate) const PROGRAM_ATTRIBUTE: &str = "program";
pub(crate) const PROGRAM_VERSION_ATTRIBUTE: &str = "programVersion";
pub(crate) const FRAME: &str = "frame";
pub(crate) const ATOM: &str = "atom";
pub(crate) const SPATIAL: &str = "spatial";
pub(crate) const CELL_SPATIAL: &str = "cell_spatial";
pub(crate) const CELL_ANGULAR: &str = "cell_angular";
pub(crate) const LABEL: &str = "label";

pub(crate) const COORDINATES: &str = "coordinates";
pub(crate) const VELOCITIES: &str = "velocities";
pub(crate) const FORCES: &str = "forces";
pub(crate) const TIME: &str = "time";
pub(crate) const CELL_LENGTHS: &str = "cell_lengths";
pub(crate) const CELL_ANGLES: &str = "cell_angles";
pub(crate) const UNITS_ATTRIBUTE: &str = "units";

pub(crate) const LENGTH_UNIT: &str = "angstrom";
pub(crate) const TIME_UNIT: &str = "picosecond";
pub(crate) const VELOCITY_UNIT: &str = "angstrom/picosecond";
pub(crate) const FORCE_UNIT: &str = "kilocalorie/mole/angstrom";
pub(crate) const KILOCALORIE_TO_KILOJOULE: f64 = 4.184;
pub(crate) const ANGLE_UNIT: &str = "degree";

pub(crate) fn validate_convention(file: &NcFile) -> Result<(), AmberNetcdfError> {
    let convention = attribute_string(file, CONVENTION_ATTRIBUTE)?;
    let version = attribute_string(file, VERSION_ATTRIBUTE)?;
    let is_amber = convention
        .split(|character: char| character == ',' || character.is_whitespace())
        .any(|token| token == CONVENTION);
    if is_amber && version == VERSION {
        Ok(())
    } else {
        Err(AmberNetcdfError::InvalidConvention)
    }
}

pub(crate) fn require_variable<'a>(
    file: &'a NcFile,
    name: &'static str,
) -> Result<&'a NcVariable, AmberNetcdfError> {
    file.variable(name)
        .map_err(|_| AmberNetcdfError::MissingVariable(name))
}

pub(crate) fn validate_units(
    variable: &NcVariable,
    name: &'static str,
    expected: &'static str,
) -> Result<(), AmberNetcdfError> {
    let declared = variable
        .attribute(UNITS_ATTRIBUTE)
        .and_then(|attribute| attribute.value.as_string());
    let units = match declared {
        Some(units) => units,
        None => "<missing>".to_string(),
    };
    if units.eq_ignore_ascii_case(expected) {
        Ok(())
    } else {
        Err(AmberNetcdfError::InvalidUnits {
            variable: name,
            units,
        })
    }
}

pub(crate) fn optional_attribute_string(file: &NcFile, name: &str) -> Option<String> {
    file.global_attribute(name)
        .ok()
        .and_then(|attribute| match &attribute.value {
            NcAttrValue::Chars(value) => Some(value.clone()),
            NcAttrValue::Strings(values) if values.len() == 1 => values.first().cloned(),
            _ => None,
        })
}

fn attribute_string(file: &NcFile, name: &str) -> Result<String, AmberNetcdfError> {
    optional_attribute_string(file, name).ok_or(AmberNetcdfError::InvalidConvention)
}
