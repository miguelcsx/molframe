//! CCD-governed chemical interactions delegated to native Rust kernels.

use crate::core::execution::{PyExecutionContext, default_context};
use crate::geometry::PyEigenOptions;
use crate::graph::PySpatialBackend;
use crate::structure::PyStructure;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "HydrogenBondOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHydrogenBondOptions(pdbiox::analysis::HydrogenBondOptions);

#[pymethods]
impl PyHydrogenBondOptions {
    #[new]
    #[pyo3(signature = (
        maximum_distance,
        minimum_angle,
        *,
        backend=PySpatialBackend::Auto,
        periodic=false,
    ))]
    fn new(
        maximum_distance: f32,
        minimum_angle: f64,
        backend: PySpatialBackend,
        periodic: bool,
    ) -> Self {
        Self(pdbiox::analysis::HydrogenBondOptions {
            maximum_donor_acceptor_distance: maximum_distance,
            minimum_angle_degrees: minimum_angle,
            backend: backend.into(),
            periodic,
        })
    }
}

impl PyHydrogenBondOptions {
    pub(crate) const fn native(&self) -> pdbiox::analysis::HydrogenBondOptions {
        self.0
    }
}

#[pyclass(name = "HydrogenBond", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyHydrogenBond {
    #[pyo3(get)]
    donor: u32,
    #[pyo3(get)]
    hydrogen: u32,
    #[pyo3(get)]
    acceptor: u32,
    #[pyo3(get)]
    donor_acceptor_distance: f32,
    #[pyo3(get)]
    hydrogen_acceptor_distance: f32,
    #[pyo3(get)]
    angle: f64,
}

#[pyclass(name = "SaltBridge", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySaltBridge {
    #[pyo3(get)]
    anion: u32,
    #[pyo3(get)]
    cation: u32,
    #[pyo3(get)]
    distance: f32,
}

#[pyclass(name = "PiStackingOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPiStackingOptions(pdbiox::analysis::PiStackingOptions);

#[pymethods]
impl PyPiStackingOptions {
    #[new]
    fn new(
        maximum_centre_distance: f32,
        maximum_parallel_angle: f64,
        minimum_perpendicular_angle: f64,
        plane_fit: &PyEigenOptions,
    ) -> Self {
        Self(pdbiox::analysis::PiStackingOptions {
            maximum_centre_distance,
            maximum_parallel_angle,
            minimum_perpendicular_angle,
            plane_fit: plane_fit.inner,
        })
    }
}

impl PyPiStackingOptions {
    pub(crate) const fn native(&self) -> pdbiox::analysis::PiStackingOptions {
        self.0
    }

    pub(crate) const fn from_native_parts(value: pdbiox::analysis::PiStackingOptions) -> Self {
        Self(value)
    }
}

#[pyclass(name = "CationPiOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCationPiOptions(pdbiox::analysis::CationPiOptions);

#[pymethods]
impl PyCationPiOptions {
    #[new]
    fn new(maximum_distance: f32, maximum_face_angle: f64, plane_fit: &PyEigenOptions) -> Self {
        Self(pdbiox::analysis::CationPiOptions {
            maximum_distance,
            maximum_face_angle,
            plane_fit: plane_fit.inner,
        })
    }
}

impl PyCationPiOptions {
    pub(crate) const fn native(&self) -> pdbiox::analysis::CationPiOptions {
        self.0
    }

    pub(crate) const fn from_native_parts(value: pdbiox::analysis::CationPiOptions) -> Self {
        Self(value)
    }
}

#[pyclass(name = "StackingKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyStackingKind {
    Parallel,
    TShaped,
}

#[pyclass(name = "PiStacking", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPiStacking {
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
    #[pyo3(get)]
    centre_distance: f32,
    #[pyo3(get)]
    angle: f64,
    #[pyo3(get)]
    kind: PyStackingKind,
}

#[pyclass(name = "CationPi", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCationPi {
    #[pyo3(get)]
    cation_residue: u32,
    #[pyo3(get)]
    ring_residue: u32,
    #[pyo3(get)]
    distance: f32,
}

#[pyclass(name = "WaterBridge", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyWaterBridge {
    #[pyo3(get)]
    water: u32,
    #[pyo3(get)]
    first: u32,
    #[pyo3(get)]
    second: u32,
}

#[pyclass(name = "WaterBridgeOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyWaterBridgeOptions(pub(crate) PyHydrogenBondOptions);

#[pymethods]
impl PyWaterBridgeOptions {
    #[new]
    fn new(hydrogen_bonds: &PyHydrogenBondOptions) -> Self {
        Self(*hydrogen_bonds)
    }

    #[getter]
    fn hydrogen_bonds(&self) -> PyHydrogenBondOptions {
        self.0
    }
}

impl PyWaterBridgeOptions {
    pub(crate) const fn native(&self) -> pdbiox::analysis::WaterBridgeOptions {
        pdbiox::analysis::WaterBridgeOptions {
            hydrogen_bonds: self.0.native(),
        }
    }

    pub(crate) const fn from_native_parts(value: pdbiox::analysis::WaterBridgeOptions) -> Self {
        Self(PyHydrogenBondOptions(value.hydrogen_bonds))
    }
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (options, context=None))]
    fn hydrogen_bonds(
        &self,
        py: Python<'_>,
        options: &PyHydrogenBondOptions,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<PyHydrogenBond>> {
        let structure = self.structure().clone();
        let options = options.0;
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || pdbiox::analysis::hydrogen_bonds(&structure, options, &context))
            .map(|values| values.into_iter().map(PyHydrogenBond::from).collect())
            .map_err(value_error)
    }

    #[pyo3(signature = (maximum_distance, *, backend=PySpatialBackend::Auto, context=None))]
    fn salt_bridges(
        &self,
        py: Python<'_>,
        maximum_distance: f32,
        backend: PySpatialBackend,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<PySaltBridge>> {
        let structure = self.structure().clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::salt_bridges(&structure, maximum_distance, backend.into(), &context)
        })
        .map(|values| values.into_iter().map(PySaltBridge::from).collect())
        .map_err(value_error)
    }

    fn pi_stacking(
        &self,
        py: Python<'_>,
        options: &PyPiStackingOptions,
    ) -> PyResult<Vec<PyPiStacking>> {
        let structure = self.structure().clone();
        let options = options.0;
        py.detach(move || pdbiox::analysis::pi_stacking(&structure, options))
            .map(|values| values.into_iter().map(PyPiStacking::from).collect())
            .map_err(value_error)
    }

    fn cation_pi(&self, py: Python<'_>, options: &PyCationPiOptions) -> PyResult<Vec<PyCationPi>> {
        let structure = self.structure().clone();
        let options = options.0;
        py.detach(move || pdbiox::analysis::cation_pi(&structure, options))
            .map(|values| values.into_iter().map(PyCationPi::from).collect())
            .map_err(value_error)
    }

    #[pyo3(signature = (options, context=None))]
    fn water_bridges(
        &self,
        py: Python<'_>,
        options: &PyHydrogenBondOptions,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<PyWaterBridge>> {
        let structure = self.structure().clone();
        let options = options.0;
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::water_bridges(
                &structure,
                pdbiox::analysis::WaterBridgeOptions {
                    hydrogen_bonds: options,
                },
                &context,
            )
        })
        .map(|values| values.into_iter().map(PyWaterBridge::from).collect())
        .map_err(value_error)
    }
}

/// Runs the native hydrogen-bond kernel without exposing a Python atom loop.
#[pyfunction]
#[pyo3(signature = (structure, options, context=None))]
pub(crate) fn hydrogen_bonds(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyHydrogenBondOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<PyHydrogenBond>> {
    let structure = structure.structure().clone();
    let options = options.0;
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || pdbiox::analysis::hydrogen_bonds(&structure, options, &context))
        .map(|values| values.into_iter().map(PyHydrogenBond::from).collect())
        .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, maximum_distance, backend, context=None))]
pub(crate) fn salt_bridges(
    py: Python<'_>,
    structure: &PyStructure,
    maximum_distance: f32,
    backend: PySpatialBackend,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<PySaltBridge>> {
    let structure = structure.structure().clone();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::salt_bridges(&structure, maximum_distance, backend.into(), &context)
    })
    .map(|values| values.into_iter().map(PySaltBridge::from).collect())
    .map_err(value_error)
}

#[pyfunction]
pub(crate) fn pi_stacking(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyPiStackingOptions,
) -> PyResult<Vec<PyPiStacking>> {
    let structure = structure.structure().clone();
    let options = options.0;
    py.detach(move || pdbiox::analysis::pi_stacking(&structure, options))
        .map(|values| values.into_iter().map(PyPiStacking::from).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn cation_pi(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyCationPiOptions,
) -> PyResult<Vec<PyCationPi>> {
    let structure = structure.structure().clone();
    let options = options.0;
    py.detach(move || pdbiox::analysis::cation_pi(&structure, options))
        .map(|values| values.into_iter().map(PyCationPi::from).collect())
        .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, options, context=None))]
pub(crate) fn water_bridges(
    py: Python<'_>,
    structure: &PyStructure,
    options: &PyWaterBridgeOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<PyWaterBridge>> {
    let structure = structure.structure().clone();
    let options = options.0.0;
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::water_bridges(
            &structure,
            pdbiox::analysis::WaterBridgeOptions {
                hydrogen_bonds: options,
            },
            &context,
        )
    })
    .map(|values| values.into_iter().map(PyWaterBridge::from).collect())
    .map_err(value_error)
}

impl From<pdbiox::analysis::HydrogenBond> for PyHydrogenBond {
    fn from(value: pdbiox::analysis::HydrogenBond) -> Self {
        Self {
            donor: value.donor.get(),
            hydrogen: value.hydrogen.get(),
            acceptor: value.acceptor.get(),
            donor_acceptor_distance: value.donor_acceptor_distance,
            hydrogen_acceptor_distance: value.hydrogen_acceptor_distance,
            angle: value.angle_degrees,
        }
    }
}

impl From<pdbiox::analysis::SaltBridge> for PySaltBridge {
    fn from(value: pdbiox::analysis::SaltBridge) -> Self {
        Self {
            anion: value.anion.get(),
            cation: value.cation.get(),
            distance: value.distance,
        }
    }
}

impl From<pdbiox::analysis::PiStacking> for PyPiStacking {
    fn from(value: pdbiox::analysis::PiStacking) -> Self {
        Self {
            first: value.first.get(),
            second: value.second.get(),
            centre_distance: value.centre_distance,
            angle: value.angle,
            kind: value.kind.into(),
        }
    }
}

impl From<pdbiox::analysis::StackingKind> for PyStackingKind {
    fn from(value: pdbiox::analysis::StackingKind) -> Self {
        match value {
            pdbiox::analysis::StackingKind::Parallel => Self::Parallel,
            pdbiox::analysis::StackingKind::TShaped => Self::TShaped,
        }
    }
}

impl From<pdbiox::analysis::CationPi> for PyCationPi {
    fn from(value: pdbiox::analysis::CationPi) -> Self {
        Self {
            cation_residue: value.cation_residue.get(),
            ring_residue: value.ring_residue.get(),
            distance: value.distance,
        }
    }
}

impl From<pdbiox::analysis::WaterBridge> for PyWaterBridge {
    fn from(value: pdbiox::analysis::WaterBridge) -> Self {
        Self {
            water: value.water.get(),
            first: value.first.get(),
            second: value.second.get(),
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
