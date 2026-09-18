//! Native-to-Python conversions for structure analyses.

use super::structural::{
    PyBasePair, PyContactMap, PyHalfSphereExposure, PyNativeContacts, PyResidueContact,
};

impl From<molframe::analysis::BasePair> for PyBasePair {
    fn from(value: molframe::analysis::BasePair) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            hydrogen_bond_count: value.hydrogen_bond_count,
            closest_distance: value.closest_distance,
        }
    }
}

impl From<molframe::analysis::ContactMap> for PyContactMap {
    fn from(value: molframe::analysis::ContactMap) -> Self {
        Self {
            residue_count: value.residue_count(),
            contacts: value
                .contacts()
                .iter()
                .copied()
                .map(PyResidueContact::from)
                .collect(),
        }
    }
}

impl From<molframe::analysis::ResidueContact> for PyResidueContact {
    fn from(value: molframe::analysis::ResidueContact) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            min_distance: value.min_distance,
        }
    }
}

impl From<molframe::analysis::HalfSphereExposure> for PyHalfSphereExposure {
    fn from(value: molframe::analysis::HalfSphereExposure) -> Self {
        Self {
            residue: value.residue.get(),
            upper: value.upper,
            lower: value.lower,
        }
    }
}

impl From<molframe::analysis::NativeContacts> for PyNativeContacts {
    fn from(value: molframe::analysis::NativeContacts) -> Self {
        Self {
            native: value.native,
            kept: value.kept,
            fraction: value.fraction,
        }
    }
}
