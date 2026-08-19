//! Native-to-Python chemistry value conversions.

use super::{PyElementProperties, PyIonicRadius, PyIonicSpin, PyRadiusSet};

impl From<PyRadiusSet> for pdbiox::RadiusSet {
    fn from(value: PyRadiusSet) -> Self {
        match value {
            PyRadiusSet::Bondi => Self::Bondi,
            PyRadiusSet::AmberUnited => Self::AmberUnited,
            PyRadiusSet::Charmm => Self::Charmm,
            PyRadiusSet::Alvarez => Self::Alvarez,
        }
    }
}

impl From<pdbiox::RadiusSet> for PyRadiusSet {
    fn from(value: pdbiox::RadiusSet) -> Self {
        match value {
            pdbiox::RadiusSet::Bondi => Self::Bondi,
            pdbiox::RadiusSet::AmberUnited => Self::AmberUnited,
            pdbiox::RadiusSet::Charmm => Self::Charmm,
            pdbiox::RadiusSet::Alvarez => Self::Alvarez,
        }
    }
}

impl From<pdbiox::ElementProperties> for PyElementProperties {
    fn from(value: pdbiox::ElementProperties) -> Self {
        Self {
            atomic_weight: value.atomic_weight,
            covalent_radius: value.covalent_radius,
            electronegativity: value.electronegativity,
            valence_electrons: value.valence_electrons,
            period: value.period,
            group: value.group,
        }
    }
}

impl From<pdbiox::IonicRadius> for PyIonicRadius {
    fn from(value: pdbiox::IonicRadius) -> Self {
        Self {
            charge: value.charge,
            coordination: value.coordination.into(),
            spin: value.spin.into(),
            ionic_radius: value.ionic_radius,
            crystal_radius: value.crystal_radius,
            most_reliable: value.most_reliable,
        }
    }
}

impl From<pdbiox::IonicSpin> for PyIonicSpin {
    fn from(value: pdbiox::IonicSpin) -> Self {
        match value {
            pdbiox::IonicSpin::Unspecified => Self::Unspecified,
            pdbiox::IonicSpin::High => Self::High,
            pdbiox::IonicSpin::Low => Self::Low,
        }
    }
}
