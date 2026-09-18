//! Native-to-Python distribution conversions.

use super::distributions::{PyLeaflet, PyRadialBin};

impl From<molframe::analysis::RadialBin> for PyRadialBin {
    fn from(bin: molframe::analysis::RadialBin) -> Self {
        Self {
            lower: bin.lower,
            upper: bin.upper,
            count: bin.count,
            distribution: bin.distribution,
        }
    }
}

impl From<molframe::analysis::Leaflet> for PyLeaflet {
    fn from(value: molframe::analysis::Leaflet) -> Self {
        Self { sites: value.sites }
    }
}
