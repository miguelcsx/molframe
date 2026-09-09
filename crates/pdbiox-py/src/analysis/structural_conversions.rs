//! Native-to-Python conversions for structure analyses.

use super::structural::{
    PyBasePair, PyContactMap, PyHalfSphereExposure, PyNativeContacts, PyResidueContact,
};

impl From<pdbiox::analysis::BasePair> for PyBasePair {
    fn from(value: pdbiox::analysis::BasePair) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            hydrogen_bond_count: value.hydrogen_bond_count,
            closest_distance: value.closest_distance,
        }
    }
}

impl From<pdbiox::analysis::ContactMap> for PyContactMap {
    fn from(value: pdbiox::analysis::ContactMap) -> Self {
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

impl From<pdbiox::analysis::ResidueContact> for PyResidueContact {
    fn from(value: pdbiox::analysis::ResidueContact) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            min_distance: value.min_distance,
        }
    }
}

impl From<pdbiox::analysis::HalfSphereExposure> for PyHalfSphereExposure {
    fn from(value: pdbiox::analysis::HalfSphereExposure) -> Self {
        Self {
            residue: value.residue.get(),
            upper: value.upper,
            lower: value.lower,
        }
    }
}

impl From<pdbiox::analysis::NativeContacts> for PyNativeContacts {
    fn from(value: pdbiox::analysis::NativeContacts) -> Self {
        Self {
            native: value.native,
            kept: value.kept,
            fraction: value.fraction,
        }
    }
}
