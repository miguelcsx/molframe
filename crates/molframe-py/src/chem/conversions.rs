//! Native-to-Python chemistry value conversions.

use super::{PyElementProperties, PyIonicRadius, PyIonicSpin, PyRadiusSet};

impl From<PyRadiusSet> for molframe::RadiusSet {
    fn from(value: PyRadiusSet) -> Self {
        match value {
            PyRadiusSet::Bondi => Self::Bondi,
            PyRadiusSet::AmberUnited => Self::AmberUnited,
            PyRadiusSet::Charmm => Self::Charmm,
            PyRadiusSet::Alvarez => Self::Alvarez,
        }
    }
}

impl From<molframe::RadiusSet> for PyRadiusSet {
    fn from(value: molframe::RadiusSet) -> Self {
        match value {
            molframe::RadiusSet::Bondi => Self::Bondi,
            molframe::RadiusSet::AmberUnited => Self::AmberUnited,
            molframe::RadiusSet::Charmm => Self::Charmm,
            molframe::RadiusSet::Alvarez => Self::Alvarez,
        }
    }
}

impl From<molframe::ElementProperties> for PyElementProperties {
    fn from(value: molframe::ElementProperties) -> Self {
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

impl From<molframe::IonicRadius> for PyIonicRadius {
    fn from(value: molframe::IonicRadius) -> Self {
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

impl From<molframe::IonicSpin> for PyIonicSpin {
    fn from(value: molframe::IonicSpin) -> Self {
        match value {
            molframe::IonicSpin::Unspecified => Self::Unspecified,
            molframe::IonicSpin::High => Self::High,
            molframe::IonicSpin::Low => Self::Low,
        }
    }
}
