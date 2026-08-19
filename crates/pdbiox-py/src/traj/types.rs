//! Typed Python controls for trajectory formats and explicit units.

use pdbiox::traj::{
    AmberNetcdfPrecision, AmberNetcdfWriteOptions, DcdEndian, DcdWriteOptions, GsdOptions,
    H5mdOptions, H5mdUnitSystem, TngCompression, TngWriteOptions, TrajectoryFormat,
    TrajectoryWriteOptions, TrrPrecision, TrrWriteOptions, XtcWriteOptions,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "TrajectoryFormat", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyTrajectoryFormat {
    Xtc,
    Trr,
    Dcd,
    AmberNetcdf,
    Tng,
    Gsd,
    H5md,
    Trz,
    Namd,
    AmberRestart,
    AmberAscii,
    Gro,
    Xyz,
    Aims,
    Txyz,
    DlPolyConfig,
    DlPolyHistory,
    CharmmCard,
    Gamess,
    LammpsDump,
    Gromos11,
    Dms,
}

impl From<PyTrajectoryFormat> for TrajectoryFormat {
    fn from(value: PyTrajectoryFormat) -> Self {
        match value {
            PyTrajectoryFormat::Xtc => Self::Xtc,
            PyTrajectoryFormat::Trr => Self::Trr,
            PyTrajectoryFormat::Dcd => Self::Dcd,
            PyTrajectoryFormat::AmberNetcdf => Self::AmberNetcdf,
            PyTrajectoryFormat::Tng => Self::Tng,
            PyTrajectoryFormat::Gsd => Self::Gsd,
            PyTrajectoryFormat::H5md => Self::H5md,
            PyTrajectoryFormat::Trz => Self::Trz,
            PyTrajectoryFormat::Namd => Self::Namd,
            PyTrajectoryFormat::AmberRestart => Self::AmberRestart,
            PyTrajectoryFormat::AmberAscii => Self::AmberAscii,
            PyTrajectoryFormat::Gro => Self::Gro,
            PyTrajectoryFormat::Xyz => Self::Xyz,
            PyTrajectoryFormat::Aims => Self::Aims,
            PyTrajectoryFormat::Txyz => Self::Txyz,
            PyTrajectoryFormat::DlPolyConfig => Self::DlPolyConfig,
            PyTrajectoryFormat::DlPolyHistory => Self::DlPolyHistory,
            PyTrajectoryFormat::CharmmCard => Self::CharmmCard,
            PyTrajectoryFormat::Gamess => Self::Gamess,
            PyTrajectoryFormat::LammpsDump => Self::LammpsDump,
            PyTrajectoryFormat::Gromos11 => Self::Gromos11,
            PyTrajectoryFormat::Dms => Self::Dms,
        }
    }
}

impl TryFrom<TrajectoryFormat> for PyTrajectoryFormat {
    type Error = PyErr;

    fn try_from(value: TrajectoryFormat) -> Result<Self, Self::Error> {
        Ok(match value {
            TrajectoryFormat::Xtc => Self::Xtc,
            TrajectoryFormat::Trr => Self::Trr,
            TrajectoryFormat::Dcd => Self::Dcd,
            TrajectoryFormat::AmberNetcdf => Self::AmberNetcdf,
            TrajectoryFormat::Tng => Self::Tng,
            TrajectoryFormat::Gsd => Self::Gsd,
            TrajectoryFormat::H5md => Self::H5md,
            TrajectoryFormat::Trz => Self::Trz,
            TrajectoryFormat::Namd => Self::Namd,
            TrajectoryFormat::AmberRestart => Self::AmberRestart,
            TrajectoryFormat::AmberAscii => Self::AmberAscii,
            TrajectoryFormat::Gro => Self::Gro,
            TrajectoryFormat::Xyz => Self::Xyz,
            TrajectoryFormat::Aims => Self::Aims,
            TrajectoryFormat::Txyz => Self::Txyz,
            TrajectoryFormat::DlPolyConfig => Self::DlPolyConfig,
            TrajectoryFormat::DlPolyHistory => Self::DlPolyHistory,
            TrajectoryFormat::CharmmCard => Self::CharmmCard,
            TrajectoryFormat::Gamess => Self::Gamess,
            TrajectoryFormat::LammpsDump => Self::LammpsDump,
            TrajectoryFormat::Gromos11 => Self::Gromos11,
            TrajectoryFormat::Dms => Self::Dms,
            _ => {
                return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                    "trajectory format is newer than this binding version",
                ));
            }
        })
    }
}

#[pyclass(name = "TrajectoryUnits", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTrajectoryUnits {
    pub(super) inner: H5mdUnitSystem,
}

#[pymethods]
impl PyTrajectoryUnits {
    #[new]
    fn new(
        length: (String, f64),
        time: (String, f64),
        velocity: (String, f64),
        force: (String, f64),
    ) -> Self {
        Self {
            inner: H5mdUnitSystem {
                length_unit: length.0,
                length_to_angstrom: length.1,
                time_unit: time.0,
                time_to_picosecond: time.1,
                velocity_unit: velocity.0,
                velocity_to_angstrom_per_picosecond: velocity.1,
                force_unit: force.0,
                force_to_kilojoule_per_mole_angstrom: force.1,
            },
        }
    }

    #[staticmethod]
    fn canonical() -> Self {
        Self {
            inner: H5mdUnitSystem::canonical(),
        }
    }

    #[getter]
    fn length_unit(&self) -> &str {
        &self.inner.length_unit
    }

    #[getter]
    const fn length_to_angstrom(&self) -> f64 {
        self.inner.length_to_angstrom
    }

    #[getter]
    fn time_unit(&self) -> &str {
        &self.inner.time_unit
    }

    #[getter]
    const fn time_to_picosecond(&self) -> f64 {
        self.inner.time_to_picosecond
    }

    #[getter]
    fn velocity_unit(&self) -> &str {
        &self.inner.velocity_unit
    }

    #[getter]
    const fn velocity_to_angstrom_per_picosecond(&self) -> f64 {
        self.inner.velocity_to_angstrom_per_picosecond
    }

    #[getter]
    fn force_unit(&self) -> &str {
        &self.inner.force_unit
    }

    #[getter]
    const fn force_to_kilojoule_per_mole_angstrom(&self) -> f64 {
        self.inner.force_to_kilojoule_per_mole_angstrom
    }
}

#[pyclass(name = "TrajectoryWriteOptions", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTrajectoryWriteOptions {
    pub(super) inner: TrajectoryWriteOptions,
}

#[pymethods]
impl PyTrajectoryWriteOptions {
    #[staticmethod]
    #[pyo3(signature = (*, precision=None))]
    fn xtc(precision: Option<f32>) -> Self {
        let mut value = XtcWriteOptions::default();
        if let Some(precision) = precision {
            value.precision = precision;
        }
        Self::new(TrajectoryFormat::Xtc, |options| options.xtc = Some(value))
    }

    #[staticmethod]
    #[pyo3(signature = (*, double_precision=false))]
    fn trr(double_precision: bool) -> Self {
        let precision = if double_precision {
            TrrPrecision::Double
        } else {
            TrrPrecision::Single
        };
        Self::new(TrajectoryFormat::Trr, |options| {
            options.trr = Some(TrrWriteOptions { precision });
        })
    }

    #[staticmethod]
    #[pyo3(signature = (*, little_endian=None, title=None, start_step=None, save_interval=None, delta_akma=None))]
    fn dcd(
        little_endian: Option<bool>,
        title: Option<String>,
        start_step: Option<i32>,
        save_interval: Option<i32>,
        delta_akma: Option<f32>,
    ) -> Self {
        let mut value = DcdWriteOptions::default();
        if let Some(little_endian) = little_endian {
            value.endian = if little_endian {
                DcdEndian::Little
            } else {
                DcdEndian::Big
            };
        }
        if let Some(title) = title {
            value.title = title.into();
        }
        if let Some(start_step) = start_step {
            value.start_step = start_step;
        }
        if let Some(save_interval) = save_interval {
            value.save_interval = save_interval;
        }
        if let Some(delta_akma) = delta_akma {
            value.delta_akma = delta_akma;
        }
        Self::new(TrajectoryFormat::Dcd, |options| options.dcd = Some(value))
    }

    #[staticmethod]
    #[pyo3(signature = (*, double_precision=false))]
    fn amber_netcdf(double_precision: bool) -> Self {
        let precision = if double_precision {
            AmberNetcdfPrecision::Double
        } else {
            AmberNetcdfPrecision::Single
        };
        Self::new(TrajectoryFormat::AmberNetcdf, |options| {
            options.amber_netcdf = Some(AmberNetcdfWriteOptions {
                precision,
                ..AmberNetcdfWriteOptions::default()
            });
        })
    }

    #[staticmethod]
    #[pyo3(signature = (*, distance_unit_exponent=None, lossy_precision=None, hashes=None))]
    fn tng(
        distance_unit_exponent: Option<i64>,
        lossy_precision: Option<f64>,
        hashes: Option<bool>,
    ) -> Self {
        let mut value = TngWriteOptions::default();
        if let Some(distance_unit_exponent) = distance_unit_exponent {
            value.distance_unit_exponent = distance_unit_exponent;
        }
        if let Some(lossy_precision) = lossy_precision {
            value.compression = TngCompression::Lossy {
                precision: lossy_precision,
            };
        }
        if let Some(hashes) = hashes {
            value.hashes = hashes;
        }
        Self::new(TrajectoryFormat::Tng, |options| options.tng = Some(value))
    }

    #[staticmethod]
    fn gsd(length_to_angstrom: f64) -> PyResult<Self> {
        let value = GsdOptions::new(length_to_angstrom)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Ok(Self::new(TrajectoryFormat::Gsd, |options| {
            options.gsd = Some(value);
        }))
    }

    #[staticmethod]
    #[pyo3(signature = (*, particle_group=None, units=None))]
    fn h5md(particle_group: Option<String>, units: Option<PyRef<'_, PyTrajectoryUnits>>) -> Self {
        let mut value = H5mdOptions::default();
        if let Some(particle_group) = particle_group {
            value.particle_group = particle_group;
        }
        if let Some(units) = units {
            value.units = units.inner.clone();
        }
        Self::new(TrajectoryFormat::H5md, |options| {
            options.h5md = Some(value);
        })
    }

    #[staticmethod]
    fn trz() -> Self {
        Self::new(TrajectoryFormat::Trz, |_| {})
    }
}

impl PyTrajectoryWriteOptions {
    fn new(format: TrajectoryFormat, configure: impl FnOnce(&mut TrajectoryWriteOptions)) -> Self {
        let mut inner = TrajectoryWriteOptions {
            format: Some(format),
            ..TrajectoryWriteOptions::default()
        };
        configure(&mut inner);
        Self { inner }
    }
}
