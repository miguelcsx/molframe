//! Explicit H5MD path and unit policy.

use super::H5mdError;

const DEFAULT_PARTICLE_GROUP: &str = "all";

/// Unit declarations and conversion factors for one H5MD file.
///
/// Each factor converts one stored value into molframe canonical units:
/// ångström, picoseconds, ångström/ps, and kJ mol⁻¹ Å⁻¹. This makes H5MD's
/// intentionally open unit vocabulary explicit rather than guessed.
#[derive(Clone, Debug, PartialEq)]
pub struct H5mdUnitSystem {
    /// Exact `unit` attribute for position and box edges.
    pub length_unit: String,
    /// Stored length multiplied by this factor gives ångström.
    pub length_to_angstrom: f64,
    /// Exact `unit` attribute for time.
    pub time_unit: String,
    /// Stored time multiplied by this factor gives picoseconds.
    pub time_to_picosecond: f64,
    /// Exact `unit` attribute for velocity.
    pub velocity_unit: String,
    /// Stored velocity multiplied by this factor gives ångström/ps.
    pub velocity_to_angstrom_per_picosecond: f64,
    /// Exact `unit` attribute for force.
    pub force_unit: String,
    /// Stored force multiplied by this factor gives kJ mol⁻¹ Å⁻¹.
    pub force_to_kilojoule_per_mole_angstrom: f64,
}

impl H5mdUnitSystem {
    /// Canonical molframe units, written without conversion.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            length_unit: "angstrom".into(),
            length_to_angstrom: 1.0,
            time_unit: "picosecond".into(),
            time_to_picosecond: 1.0,
            velocity_unit: "angstrom/picosecond".into(),
            velocity_to_angstrom_per_picosecond: 1.0,
            force_unit: "kilojoule/mole/angstrom".into(),
            force_to_kilojoule_per_mole_angstrom: 1.0,
        }
    }

    pub(crate) fn validate(&self) -> Result<(), H5mdError> {
        let named = [
            &self.length_unit,
            &self.time_unit,
            &self.velocity_unit,
            &self.force_unit,
        ];
        let scales = [
            self.length_to_angstrom,
            self.time_to_picosecond,
            self.velocity_to_angstrom_per_picosecond,
            self.force_to_kilojoule_per_mole_angstrom,
        ];
        if named
            .iter()
            .all(|name| !name.is_empty() && !name.contains('\0'))
            && scales.iter().all(|scale| scale.is_finite() && *scale > 0.0)
        {
            Ok(())
        } else {
            Err(H5mdError::InvalidValue)
        }
    }
}

impl Default for H5mdUnitSystem {
    fn default() -> Self {
        Self::canonical()
    }
}

/// H5MD particle group and unit policy.
#[derive(Clone, Debug, PartialEq)]
pub struct H5mdOptions {
    /// Name below `/particles` to read or write.
    pub particle_group: String,
    /// Exact unit declarations and conversion factors.
    pub units: H5mdUnitSystem,
}

impl H5mdOptions {
    pub(crate) fn validate(&self) -> Result<(), H5mdError> {
        if self.particle_group.is_empty()
            || self.particle_group.contains('/')
            || self.particle_group.contains('\0')
        {
            return Err(H5mdError::InvalidValue);
        }
        self.units.validate()
    }
}

impl Default for H5mdOptions {
    fn default() -> Self {
        Self {
            particle_group: DEFAULT_PARTICLE_GROUP.into(),
            units: H5mdUnitSystem::canonical(),
        }
    }
}
