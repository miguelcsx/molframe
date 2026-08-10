//! Unified trajectory I/O errors.

/// Path dispatch, option or format-specific failure.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum TrajectoryIoError {
    /// The path has no recognised trajectory suffix and no format was selected.
    #[error("trajectory format cannot be inferred from the path")]
    UnknownFormat,
    /// A unit-agnostic format needs an explicit conversion.
    #[error("trajectory format requires explicit unit options")]
    MissingUnitOptions,
    /// A format without an atom-count header needs its paired topology size.
    #[error("trajectory format requires explicit atom-count and box options")]
    MissingTopologyOptions,
    /// Per-frame metadata is absent, inconsistent or not representable.
    #[error("trajectory metadata is inconsistent with its frames")]
    InvalidMetadata,
    /// A numeric value cannot be represented by the canonical trajectory model.
    #[error("trajectory value is outside the canonical numeric representation")]
    UnrepresentableValue,
    /// The selected format is intentionally read-only in the public contract.
    #[error("selected trajectory format is read-only")]
    ReadOnlyFormat,
    /// Filesystem access failed.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// XTC decoding or encoding failed.
    #[error(transparent)]
    Xtc(#[from] crate::XtcError),
    /// TRR decoding or encoding failed.
    #[error(transparent)]
    Trr(#[from] crate::TrrError),
    /// DCD decoding or encoding failed.
    #[error(transparent)]
    Dcd(#[from] crate::DcdError),
    /// AMBER `NetCDF` decoding or encoding failed.
    #[error(transparent)]
    AmberNetcdf(#[from] crate::AmberNetcdfError),
    /// TNG decoding or encoding failed.
    #[error(transparent)]
    Tng(#[from] crate::TngError),
    /// GSD decoding or encoding failed.
    #[error(transparent)]
    Gsd(#[from] crate::GsdError),
    /// H5MD decoding or encoding failed.
    #[error(transparent)]
    H5md(#[from] crate::H5mdError),
    /// TRZ decoding or encoding failed.
    #[error(transparent)]
    Trz(#[from] crate::TrzError),
    /// NAMD binary snapshot decoding or encoding failed.
    #[error(transparent)]
    Namd(#[from] crate::NamdError),
    /// Formatted AMBER restart decoding or encoding failed.
    #[error(transparent)]
    Amber(#[from] crate::AmberError),
    /// GRO text decoding or encoding failed.
    #[error(transparent)]
    Gro(#[from] crate::GroError),
    /// FHI-aims geometry decoding failed.
    #[error(transparent)]
    Aims(#[from] crate::AimsError),
    /// Tinker XYZ/ARC decoding or encoding failed.
    #[error(transparent)]
    Txyz(#[from] crate::TxyzError),
    /// `DL_POLY` decoding or encoding failed.
    #[error(transparent)]
    DlPoly(#[from] crate::DlPolyError),
    /// CHARMM coordinate CARD decoding or encoding failed.
    #[error(transparent)]
    Charmm(#[from] crate::CharmmError),
    /// GAMESS output decoding failed.
    #[error(transparent)]
    Gamess(#[from] crate::GamessError),
    /// LAMMPS dump decoding failed.
    #[error(transparent)]
    Lammps(#[from] crate::LammpsError),
    /// GROMOS11 decoding failed.
    #[error(transparent)]
    Gromos(#[from] crate::GromosError),
    /// DESRES DMS database decoding failed.
    #[error(transparent)]
    Dms(#[from] crate::DmsError),
}
