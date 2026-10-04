//! The words a trajectory container is spelled with.

use super::model::TrajectoryFormat;
use molframe_core::contract::{PolicyParseError, canonical_spelling};
use std::fmt;
use std::str::FromStr;

/// Every container, with its canonical spelling.
const FORMATS: &[(TrajectoryFormat, &str)] = &[
    (TrajectoryFormat::Xtc, "xtc"),
    (TrajectoryFormat::Trr, "trr"),
    (TrajectoryFormat::Dcd, "dcd"),
    (TrajectoryFormat::AmberNetcdf, "amber-netcdf"),
    (TrajectoryFormat::Tng, "tng"),
    (TrajectoryFormat::Gsd, "gsd"),
    (TrajectoryFormat::H5md, "h5md"),
    (TrajectoryFormat::Trz, "trz"),
    (TrajectoryFormat::Namd, "namd"),
    (TrajectoryFormat::AmberRestart, "amber-restart"),
    (TrajectoryFormat::AmberAscii, "amber-ascii"),
    (TrajectoryFormat::Gro, "gro"),
    (TrajectoryFormat::Xyz, "xyz"),
    (TrajectoryFormat::Aims, "aims"),
    (TrajectoryFormat::Txyz, "txyz"),
    (TrajectoryFormat::DlPolyConfig, "dl-poly-config"),
    (TrajectoryFormat::DlPolyHistory, "dl-poly-history"),
    (TrajectoryFormat::CharmmCard, "charmm-card"),
    (TrajectoryFormat::Gamess, "gamess"),
    (TrajectoryFormat::LammpsDump, "lammps-dump"),
    (TrajectoryFormat::Gromos11, "gromos11"),
    (TrajectoryFormat::Dms, "dms"),
];

impl TrajectoryFormat {
    /// The canonical spelling of this container, `"other"` for one this build
    /// does not name.
    #[must_use]
    pub fn name(self) -> &'static str {
        match FORMATS.iter().find(|(format, _)| *format == self) {
            Some((_, name)) => name,
            None => "other",
        }
    }

    /// Every container this build names, in canonical spelling.
    pub fn names() -> impl Iterator<Item = &'static str> {
        FORMATS.iter().map(|(_, name)| *name)
    }
}

impl fmt::Display for TrajectoryFormat {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

impl FromStr for TrajectoryFormat {
    type Err = PolicyParseError;

    /// An underscore is read as a hyphen; `netcdf` is accepted for
    /// `amber-netcdf`.
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let spelled = canonical_spelling(value);
        let spelled = if spelled == "netcdf" {
            "amber-netcdf"
        } else {
            spelled.as_str()
        };
        FORMATS
            .iter()
            .find_map(|(format, name)| (*name == spelled).then_some(*format))
            .ok_or_else(|| {
                PolicyParseError::new(
                    "format",
                    value,
                    &Self::names().collect::<Vec<_>>().join(", "),
                )
            })
    }
}
