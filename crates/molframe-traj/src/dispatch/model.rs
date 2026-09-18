//! Public types for declarative trajectory dispatch.

use std::path::Path;

use crate::{
    AmberNetcdfWriteOptions, DcdHeader, DcdWriteOptions, GsdOptions, H5mdOptions, Timestep,
    TngWriteOptions, TrrPrecision, TrrWriteOptions, XtcWriteOptions,
};

/// Supported coordinate-trajectory containers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum TrajectoryFormat {
    /// GROMACS compressed coordinates.
    Xtc,
    /// GROMACS full-precision trajectory.
    Trr,
    /// CHARMM/NAMD DCD.
    Dcd,
    /// AMBER `NetCDF` trajectory convention.
    AmberNetcdf,
    /// Trajectory Next Generation.
    Tng,
    /// HOOMD GSD.
    Gsd,
    /// H5MD.
    H5md,
    /// IBIsCO/YASP TRZ.
    Trz,
    /// NAMD binary coordinate snapshot.
    Namd,
    /// Formatted AMBER restart/inpcrd snapshot.
    AmberRestart,
    /// AMBER ASCII trajectory with topology supplied by the caller.
    AmberAscii,
    /// GROMACS text coordinate trajectory.
    Gro,
    /// Multi-frame XYZ coordinates with elements and comments.
    Xyz,
    /// FHI-aims `geometry.in` coordinates.
    Aims,
    /// Tinker XYZ/TXYZ or multi-frame ARC.
    Txyz,
    /// `DL_POLY` CONFIG snapshot.
    DlPolyConfig,
    /// `DL_POLY` HISTORY trajectory.
    DlPolyHistory,
    /// CHARMM coordinate CARD.
    CharmmCard,
    /// GAMESS optimization or surface output.
    Gamess,
    /// LAMMPS native text dump.
    LammpsDump,
    /// GROMOS11 block trajectory.
    Gromos11,
    /// DESRES molecular structure database.
    Dms,
}

impl TrajectoryFormat {
    /// Infers a trajectory container from a case-insensitive file suffix.
    #[must_use]
    pub fn infer(path: &Path) -> Option<Self> {
        let extension = path.extension()?.to_str()?.to_ascii_lowercase();
        match extension.as_str() {
            "xtc" => Some(Self::Xtc),
            "trr" => Some(Self::Trr),
            "dcd" => Some(Self::Dcd),
            "nc" | "ncdf" | "netcdf" => Some(Self::AmberNetcdf),
            "tng" => Some(Self::Tng),
            "gsd" => Some(Self::Gsd),
            "h5md" => Some(Self::H5md),
            "trz" => Some(Self::Trz),
            "coor" | "namdbin" => Some(Self::Namd),
            "rst7" | "inpcrd" | "restrt" => Some(Self::AmberRestart),
            "trj" | "mdcrd" | "crdbox" => Some(Self::AmberAscii),
            "gro" => Some(Self::Gro),
            "xyz" => Some(Self::Xyz),
            "aims" | "in" => Some(Self::Aims),
            "txyz" | "arc" => Some(Self::Txyz),
            "config" => Some(Self::DlPolyConfig),
            "history" => Some(Self::DlPolyHistory),
            "crd" => Some(Self::CharmmCard),
            "gms" | "log" | "out" => Some(Self::Gamess),
            "lammpsdump" => Some(Self::LammpsDump),
            "trc" => Some(Self::Gromos11),
            "dms" => Some(Self::Dms),
            _ => None,
        }
    }
}

/// Resource ceiling and format selection for pull-based trajectory readers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TrajectoryReaderOptions {
    /// Override suffix-based format inference.
    pub format: Option<TrajectoryFormat>,
    /// Maximum bytes used by reader workspaces, indexes and one output frame.
    pub memory_limit_bytes: usize,
}

impl TrajectoryReaderOptions {
    /// Conservative default suitable for long-running analysis processes.
    pub const DEFAULT_MEMORY_LIMIT_BYTES: usize = 100_000_000;
}

impl Default for TrajectoryReaderOptions {
    fn default() -> Self {
        Self {
            format: None,
            memory_limit_bytes: Self::DEFAULT_MEMORY_LIMIT_BYTES,
        }
    }
}

/// Source-specific metadata required for faithful inspection and rewriting.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum FormatMetadata {
    /// No additional container controls.
    None,
    /// Per-frame inverse-nanometre coordinate quantization.
    Xtc {
        /// Inverse-nanometre quantization for each frame.
        precision: Vec<f32>,
    },
    /// Per-frame numeric representation.
    Trr {
        /// Floating-point representation for each frame.
        precision: Vec<TrrPrecision>,
    },
    /// Complete DCD header.
    Dcd(DcdHeader),
    /// AMBER numeric representation and producer attributes.
    AmberNetcdf(crate::AmberNetcdfMetadata),
    /// TNG unit and compression header values.
    Tng {
        /// Base-ten exponent of the source length unit in metres.
        distance_unit_exponent: i64,
        /// Header coordinate precision.
        compression_precision: f64,
        /// Position codec retained from the source file.
        compression: crate::TngCompression,
    },
    /// Explicit conversion used for unit-agnostic GSD coordinates.
    Gsd(GsdOptions),
    /// H5MD particle group, units and creator identity.
    H5md(crate::H5mdMetadata),
    /// TRZ title and force-stream presence.
    Trz {
        /// Fixed-width title without padding.
        title: Box<str>,
        /// Whether frames carry forces.
        has_forces: bool,
    },
    /// NAMD byte order retained from its binary snapshot.
    Namd(crate::NamdEndian),
    /// AMBER title and explicit scalar layout.
    AmberRestart {
        /// Title line.
        title: Box<str>,
        /// Coordinate/velocity/cell layout.
        layout: crate::AmberRestartLayout,
    },
    /// AMBER ASCII title and explicit topology interpretation.
    AmberAscii {
        /// Title line retained verbatim.
        title: Box<str>,
        /// Atom count supplied by the topology.
        atom_count: usize,
        /// Whether every frame carries three orthorhombic box lengths.
        periodic_box: bool,
    },
    /// Complete GRO records carrying topology, titles and periodic boxes.
    Gro(Vec<crate::GroFrame>),
    /// Complete XYZ records carrying elements and comments.
    Xyz(Vec<crate::XyzFrame>),
    /// Complete FHI-aims species and lattice metadata.
    Aims(crate::AimsGeometry),
    /// Complete Tinker atom types and connectivity.
    Txyz(Vec<crate::TxyzFrame>),
    /// Complete `DL_POLY` CONFIG metadata.
    DlPolyConfig(crate::DlPolyConfig),
    /// Complete `DL_POLY` HISTORY metadata.
    DlPolyHistory(crate::DlPolyHistory),
    /// Complete CHARMM atom-card metadata.
    CharmmCard(crate::CharmmCard),
    /// GAMESS atoms, run type and scalar frame metadata.
    Gamess(crate::GamessTrajectory),
    /// GROMOS title and boundary convention.
    Gromos11(crate::GromosTrajectory),
    /// Complete DESRES topology, annotations, coordinates and cell.
    Dms(Box<crate::DmsSystem>),
}

/// Cross-format metadata that is not part of [`Timestep`].
#[derive(Clone, Debug, PartialEq)]
pub struct TrajectoryMetadata {
    /// Simulation steps in frame order, when the format carries them.
    pub steps: Option<Vec<i64>>,
    /// Container-specific controls.
    pub format: FormatMetadata,
}

impl Default for TrajectoryMetadata {
    fn default() -> Self {
        Self {
            steps: None,
            format: FormatMetadata::None,
        }
    }
}

/// Normalized trajectory with faithful source metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct TrajectoryData {
    /// Detected or selected container.
    pub format: TrajectoryFormat,
    /// Frames in source order and canonical molframe units.
    pub frames: Vec<Timestep>,
    /// Step and container metadata.
    pub metadata: TrajectoryMetadata,
}

impl TrajectoryData {
    /// Selects frames in caller order while keeping per-frame metadata aligned.
    ///
    /// Indices may repeat. This makes extraction, reordering and duplication one
    /// declarative operation without format-specific metadata handling in
    /// callers.
    ///
    /// # Errors
    ///
    /// Returns [`crate::TrajectoryIoError::InvalidMetadata`] when an index is
    /// outside the trajectory or existing per-frame metadata is misaligned.
    pub fn select_frames(&self, indices: &[usize]) -> Result<Self, crate::TrajectoryIoError> {
        let frame_count = self.frames.len();
        if self
            .metadata
            .steps
            .as_ref()
            .is_some_and(|steps| steps.len() != frame_count)
            || !metadata_is_aligned(&self.metadata.format, frame_count)
        {
            return Err(crate::TrajectoryIoError::InvalidMetadata);
        }
        let frames = select(&self.frames, indices)?;
        let steps = self
            .metadata
            .steps
            .as_ref()
            .map(|steps| select(steps, indices))
            .transpose()?;
        let format = match &self.metadata.format {
            FormatMetadata::Xtc { precision } => FormatMetadata::Xtc {
                precision: select(precision, indices)?,
            },
            FormatMetadata::Trr { precision } => FormatMetadata::Trr {
                precision: select(precision, indices)?,
            },
            FormatMetadata::Gro(records) => FormatMetadata::Gro(select(records, indices)?),
            FormatMetadata::Xyz(records) => FormatMetadata::Xyz(select(records, indices)?),
            FormatMetadata::Txyz(records) => FormatMetadata::Txyz(select(records, indices)?),
            FormatMetadata::DlPolyHistory(source) => {
                let mut selected = source.clone();
                selected.frames = select(&source.frames, indices)?;
                FormatMetadata::DlPolyHistory(selected)
            }
            FormatMetadata::Gamess(source) => {
                let mut selected = source.clone();
                selected.frames = select(&source.frames, indices)?;
                FormatMetadata::Gamess(selected)
            }
            FormatMetadata::Gromos11(source) => {
                let mut selected = source.clone();
                selected.frames = select(&source.frames, indices)?;
                selected.boundaries = select(&source.boundaries, indices)?;
                FormatMetadata::Gromos11(selected)
            }
            other => other.clone(),
        };
        Ok(Self {
            format: self.format,
            frames,
            metadata: TrajectoryMetadata { steps, format },
        })
    }
}

fn metadata_is_aligned(metadata: &FormatMetadata, frame_count: usize) -> bool {
    match metadata {
        FormatMetadata::Xtc { precision } => precision.len() == frame_count,
        FormatMetadata::Trr { precision } => precision.len() == frame_count,
        FormatMetadata::Gro(records) => records.len() == frame_count,
        FormatMetadata::Xyz(records) => records.len() == frame_count,
        FormatMetadata::Txyz(records) => records.len() == frame_count,
        FormatMetadata::DlPolyHistory(source) => source.frames.len() == frame_count,
        FormatMetadata::Gamess(source) => source.frames.len() == frame_count,
        FormatMetadata::Gromos11(source) => {
            source.frames.len() == frame_count && source.boundaries.len() == frame_count
        }
        _ => true,
    }
}

fn select<T: Clone>(values: &[T], indices: &[usize]) -> Result<Vec<T>, crate::TrajectoryIoError> {
    indices
        .iter()
        .map(|index| {
            values
                .get(*index)
                .cloned()
                .ok_or(crate::TrajectoryIoError::InvalidMetadata)
        })
        .collect()
}

/// Explicit controls for path-based reading.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TrajectoryReadOptions {
    /// Override suffix-based format inference.
    pub format: Option<TrajectoryFormat>,
    /// Required conversion for unit-agnostic GSD files.
    pub gsd: Option<GsdOptions>,
    /// H5MD particle group and exact units.
    pub h5md: H5mdOptions,
    /// Explicit AMBER restart layout; `Auto` refuses ambiguous tails.
    pub amber_restart_layout: crate::AmberRestartLayout,
    /// Required atom count and box convention for headerless AMBER ASCII frames.
    pub amber_ascii: Option<AmberAsciiReadOptions>,
}

/// Explicit topology information required to decode AMBER ASCII trajectories.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AmberAsciiReadOptions {
    /// Number of atoms in the paired topology.
    pub atom_count: usize,
    /// Whether each frame ends with three orthorhombic box lengths.
    pub periodic_box: bool,
}

/// Explicit controls required when creating a TRZ container.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TrzWriteOptions {
    /// File title, limited by TRZ to 80 bytes.
    pub title: Box<str>,
}

/// Explicit controls for path-based writing.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TrajectoryWriteOptions {
    /// Override suffix-based format inference.
    pub format: Option<TrajectoryFormat>,
    /// XTC quantization override.
    pub xtc: Option<XtcWriteOptions>,
    /// TRR numeric representation override.
    pub trr: Option<TrrWriteOptions>,
    /// DCD header override.
    pub dcd: Option<DcdWriteOptions>,
    /// AMBER `NetCDF` precision override.
    pub amber_netcdf: Option<AmberNetcdfWriteOptions>,
    /// TNG unit/compression override.
    pub tng: Option<TngWriteOptions>,
    /// GSD explicit unit conversion.
    pub gsd: Option<GsdOptions>,
    /// H5MD particle group and exact units.
    pub h5md: Option<H5mdOptions>,
    /// TRZ container metadata; required when the source is not TRZ.
    pub trz: Option<TrzWriteOptions>,
    /// Explicit byte order for a new NAMD snapshot.
    pub namd: Option<crate::NamdEndian>,
}
