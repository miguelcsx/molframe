//! Python adaptation of the native IMD v2 TCP client.

use crate::errors::ImdError as PyImdError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::net::TcpStream;
use std::time::Duration;

#[pyclass(name = "ImdLimits", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyImdLimits {
    #[pyo3(get)]
    pub(crate) max_atoms: usize,
}

#[pymethods]
impl PyImdLimits {
    #[new]
    #[pyo3(signature = (max_atoms=None))]
    fn new(max_atoms: Option<usize>) -> PyResult<Self> {
        let max_atoms = match max_atoms {
            Some(value) => value,
            None => {
                pdbiox::traj::ImdConnectionOptions::default()
                    .limits
                    .max_atoms
            }
        };
        if max_atoms == 0 {
            return Err(PyValueError::new_err("IMD max_atoms must be positive"));
        }
        Ok(Self { max_atoms })
    }
}

impl From<PyImdLimits> for pdbiox::traj::ImdLimits {
    fn from(value: PyImdLimits) -> Self {
        Self {
            max_atoms: value.max_atoms,
        }
    }
}

impl From<pdbiox::traj::ImdLimits> for PyImdLimits {
    fn from(value: pdbiox::traj::ImdLimits) -> Self {
        Self {
            max_atoms: value.max_atoms,
        }
    }
}

#[pyclass(name = "ImdConnectionOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyImdConnectionOptions {
    inner: pdbiox::traj::ImdConnectionOptions,
}

#[pymethods]
impl PyImdConnectionOptions {
    #[new]
    #[pyo3(signature = (limits=None, *, read_timeout=None, write_timeout=None, no_delay=true))]
    fn new(
        limits: Option<PyImdLimits>,
        read_timeout: Option<f64>,
        write_timeout: Option<f64>,
        no_delay: bool,
    ) -> PyResult<Self> {
        let defaults = pdbiox::traj::ImdConnectionOptions::default();
        Ok(Self {
            inner: pdbiox::traj::ImdConnectionOptions {
                limits: limits.map_or(defaults.limits, Into::into),
                read_timeout: timeout(read_timeout)?,
                write_timeout: timeout(write_timeout)?,
                no_delay,
            },
        })
    }

    #[getter]
    fn limits(&self) -> PyImdLimits {
        self.inner.limits.into()
    }

    #[getter]
    fn read_timeout(&self) -> Option<f64> {
        self.inner.read_timeout.map(|value| value.as_secs_f64())
    }

    #[getter]
    fn write_timeout(&self) -> Option<f64> {
        self.inner.write_timeout.map(|value| value.as_secs_f64())
    }

    #[getter]
    const fn no_delay(&self) -> bool {
        self.inner.no_delay
    }
}

fn timeout(seconds: Option<f64>) -> PyResult<Option<Duration>> {
    seconds
        .map(|value| {
            if !value.is_finite() || value < 0.0 {
                return Err(PyValueError::new_err(
                    "IMD timeouts must be finite non-negative seconds",
                ));
            }
            Duration::try_from_secs_f64(value)
                .map(Some)
                .map_err(|_| PyValueError::new_err("IMD timeout is outside Duration range"))
        })
        .transpose()
        .map(Option::flatten)
}

#[pyclass(name = "ImdPeerEndian", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyImdPeerEndian {
    Little,
    Big,
}

impl From<pdbiox::traj::ImdPeerEndian> for PyImdPeerEndian {
    fn from(value: pdbiox::traj::ImdPeerEndian) -> Self {
        match value {
            pdbiox::traj::ImdPeerEndian::Little => Self::Little,
            pdbiox::traj::ImdPeerEndian::Big => Self::Big,
        }
    }
}

#[pyclass(name = "ImdEnergies", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyImdEnergies {
    #[pyo3(get)]
    step: i32,
    #[pyo3(get)]
    temperature: f32,
    #[pyo3(get)]
    total: f32,
    #[pyo3(get)]
    potential: f32,
    #[pyo3(get)]
    van_der_waals: f32,
    #[pyo3(get)]
    electrostatic: f32,
    #[pyo3(get)]
    bond: f32,
    #[pyo3(get)]
    angle: f32,
    #[pyo3(get)]
    dihedral: f32,
    #[pyo3(get)]
    improper: f32,
}

#[pymethods]
impl PyImdEnergies {
    #[new]
    #[pyo3(signature = (step, temperature, total, potential, van_der_waals, electrostatic, bond, angle, dihedral, improper))]
    fn new(
        step: i32,
        temperature: f32,
        total: f32,
        potential: f32,
        van_der_waals: f32,
        electrostatic: f32,
        bond: f32,
        angle: f32,
        dihedral: f32,
        improper: f32,
    ) -> Self {
        Self {
            step,
            temperature,
            total,
            potential,
            van_der_waals,
            electrostatic,
            bond,
            angle,
            dihedral,
            improper,
        }
    }
}

impl From<pdbiox::traj::ImdEnergies> for PyImdEnergies {
    fn from(value: pdbiox::traj::ImdEnergies) -> Self {
        Self {
            step: value.step,
            temperature: value.temperature,
            total: value.total,
            potential: value.potential,
            van_der_waals: value.van_der_waals,
            electrostatic: value.electrostatic,
            bond: value.bond,
            angle: value.angle,
            dihedral: value.dihedral,
            improper: value.improper,
        }
    }
}

impl From<PyImdEnergies> for pdbiox::traj::ImdEnergies {
    fn from(value: PyImdEnergies) -> Self {
        Self {
            step: value.step,
            temperature: value.temperature,
            total: value.total,
            potential: value.potential,
            van_der_waals: value.van_der_waals,
            electrostatic: value.electrostatic,
            bond: value.bond,
            angle: value.angle,
            dihedral: value.dihedral,
            improper: value.improper,
        }
    }
}

#[pyclass(name = "ImdForce", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyImdForce {
    #[pyo3(get)]
    atom: u32,
    #[pyo3(get)]
    force: [f32; 3],
}

#[pymethods]
impl PyImdForce {
    #[new]
    fn new(atom: u32, force: [f32; 3]) -> Self {
        Self { atom, force }
    }
}

impl From<PyImdForce> for pdbiox::traj::ImdForce {
    fn from(value: PyImdForce) -> Self {
        Self {
            atom: value.atom,
            force: value.force,
        }
    }
}

impl From<pdbiox::traj::ImdForce> for PyImdForce {
    fn from(value: pdbiox::traj::ImdForce) -> Self {
        Self {
            atom: value.atom,
            force: value.force,
        }
    }
}

#[pyclass(name = "ImdMessage", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyImdMessage {
    #[pyo3(get)]
    kind: String,
    #[pyo3(get)]
    energies: Option<PyImdEnergies>,
    #[pyo3(get)]
    coordinates: Option<Vec<[f32; 3]>>,
    #[pyo3(get)]
    forces: Option<Vec<PyImdForce>>,
    #[pyo3(get)]
    rate: Option<u32>,
}

#[pymethods]
impl PyImdMessage {
    #[staticmethod]
    fn from_disconnect() -> Self {
        Self::empty("disconnect")
    }

    #[staticmethod]
    fn from_energies(value: PyImdEnergies) -> Self {
        Self {
            kind: "energies".into(),
            energies: Some(value),
            coordinates: None,
            forces: None,
            rate: None,
        }
    }

    #[staticmethod]
    fn from_coordinates(value: Vec<[f32; 3]>) -> Self {
        Self {
            kind: "coordinates".into(),
            energies: None,
            coordinates: Some(value),
            forces: None,
            rate: None,
        }
    }

    #[staticmethod]
    fn from_go() -> Self {
        Self::empty("go")
    }

    #[staticmethod]
    fn from_kill() -> Self {
        Self::empty("kill")
    }

    #[staticmethod]
    fn from_forces(value: Vec<PyImdForce>) -> Self {
        Self {
            kind: "forces".into(),
            energies: None,
            coordinates: None,
            forces: Some(value),
            rate: None,
        }
    }

    #[staticmethod]
    fn from_pause() -> Self {
        Self::empty("pause")
    }

    #[staticmethod]
    fn from_transmission_rate(rate: u32) -> Self {
        Self {
            kind: "transmission_rate".into(),
            energies: None,
            coordinates: None,
            forces: None,
            rate: Some(rate),
        }
    }

    fn __repr__(&self) -> String {
        format!("ImdMessage(kind={:?})", self.kind)
    }
}

impl PyImdMessage {
    fn empty(kind: &str) -> Self {
        Self {
            kind: kind.into(),
            energies: None,
            coordinates: None,
            forces: None,
            rate: None,
        }
    }
}

impl From<pdbiox::traj::ImdMessage> for PyImdMessage {
    fn from(value: pdbiox::traj::ImdMessage) -> Self {
        match value {
            pdbiox::traj::ImdMessage::Disconnect => Self::from_disconnect(),
            pdbiox::traj::ImdMessage::Energies(value) => Self::from_energies(value.into()),
            pdbiox::traj::ImdMessage::Coordinates(value) => Self::from_coordinates(value),
            pdbiox::traj::ImdMessage::Go => Self::from_go(),
            pdbiox::traj::ImdMessage::Kill => Self::from_kill(),
            pdbiox::traj::ImdMessage::Forces(value) => {
                Self::from_forces(value.into_iter().map(Into::into).collect())
            }
            pdbiox::traj::ImdMessage::Pause => Self::from_pause(),
            pdbiox::traj::ImdMessage::TransmissionRate(value) => {
                Self::from_transmission_rate(value)
            }
            _ => Self::empty("unknown"),
        }
    }
}

#[pyclass(name = "ImdClient", skip_from_py_object)]
pub(crate) struct PyImdClient {
    inner: pdbiox::traj::ImdClient<TcpStream>,
}

#[pymethods]
impl PyImdClient {
    #[staticmethod]
    #[pyo3(signature = (address, *, options=None))]
    fn connect(
        py: Python<'_>,
        address: String,
        options: Option<PyImdConnectionOptions>,
    ) -> PyResult<Self> {
        let options = options.map_or_else(pdbiox::traj::ImdConnectionOptions::default, |value| {
            value.inner
        });
        py.detach(move || pdbiox::traj::ImdClient::connect_with_options(address, options))
            .map(|inner| Self { inner })
            .map_err(imd_error)
    }

    fn receive(&mut self, py: Python<'_>) -> PyResult<PyImdMessage> {
        py.detach(|| self.inner.receive())
            .map(Into::into)
            .map_err(imd_error)
    }

    fn send_forces(&mut self, py: Python<'_>, forces: Vec<PyImdForce>) -> PyResult<()> {
        let forces = forces.into_iter().map(Into::into).collect::<Vec<_>>();
        py.detach(|| self.inner.send_forces(&forces))
            .map_err(imd_error)
    }

    fn toggle_pause(&mut self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.inner.toggle_pause()).map_err(imd_error)
    }

    fn resume(&mut self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.inner.resume()).map_err(imd_error)
    }

    fn set_transmission_rate(&mut self, py: Python<'_>, steps: u32) -> PyResult<()> {
        py.detach(|| self.inner.set_transmission_rate(steps))
            .map_err(imd_error)
    }

    fn disconnect(&mut self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.inner.disconnect()).map_err(imd_error)
    }

    fn kill(&mut self, py: Python<'_>) -> PyResult<()> {
        py.detach(|| self.inner.kill()).map_err(imd_error)
    }

    #[getter]
    fn peer_endian(&self) -> PyImdPeerEndian {
        self.inner.peer_endian().into()
    }

    #[getter]
    fn closed(&self) -> bool {
        self.inner.is_closed()
    }
}

fn imd_error(error: pdbiox::traj::ImdError) -> PyErr {
    PyImdError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyImdLimits>()?;
    module.add_class::<PyImdConnectionOptions>()?;
    module.add_class::<PyImdPeerEndian>()?;
    module.add_class::<PyImdEnergies>()?;
    module.add_class::<PyImdForce>()?;
    module.add_class::<PyImdMessage>()?;
    module.add_class::<PyImdClient>()?;
    Ok(())
}
