//! Native-to-Python distribution conversions.

use super::distributions::{PyLeaflet, PyRadialBin};

impl From<pdbiox::analysis::RadialBin> for PyRadialBin {
    fn from(bin: pdbiox::analysis::RadialBin) -> Self {
        Self {
            lower: bin.lower,
            upper: bin.upper,
            count: bin.count,
            distribution: bin.distribution,
        }
    }
}

impl From<pdbiox::analysis::Leaflet> for PyLeaflet {
    fn from(value: pdbiox::analysis::Leaflet) -> Self {
        Self { sites: value.sites }
    }
}
