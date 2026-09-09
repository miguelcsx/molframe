//! Structure-level interaction analyses delegated to native Rust kernels.

use crate::chemistry::PyComponentDictionary;
use crate::contract::PyAnalysis;
use crate::core::execution::{PyExecutionContext, default_context};
use crate::graph::PySpatialBackend;
use crate::structure::PyStructure;
use numpy::PyReadonlyArray1;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "BasePairOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBasePairOptions(pdbiox::analysis::BasePairOptions);

#[pymethods]
impl PyBasePairOptions {
    #[new]
    fn new(
        maximum_donor_acceptor_distance: f32,
        minimum_angle_degrees: f64,
        minimum_hydrogen_bonds: usize,
        backend: PySpatialBackend,
        periodic: bool,
    ) -> Self {
        Self(pdbiox::analysis::BasePairOptions {
            hydrogen_bonds: pdbiox::analysis::HydrogenBondOptions {
                maximum_donor_acceptor_distance,
                minimum_angle_degrees,
                backend: backend.into(),
                periodic,
            },
            minimum_hydrogen_bonds,
        })
    }
}

impl PyBasePairOptions {
    pub(crate) const fn native(&self) -> pdbiox::analysis::BasePairOptions {
        self.0
    }

    pub(crate) const fn from_native(value: pdbiox::analysis::BasePairOptions) -> Self {
        Self(value)
    }
}

#[pyclass(name = "BasePair", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyBasePair {
    #[pyo3(get)]
    pub(super) first: u32,
    #[pyo3(get)]
    pub(super) second: u32,
    #[pyo3(get)]
    pub(super) hydrogen_bond_count: usize,
    #[pyo3(get)]
    pub(super) closest_distance: f32,
}

#[pyclass(name = "ResidueContact", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyResidueContact {
    #[pyo3(get)]
    pub(super) first: u32,
    #[pyo3(get)]
    pub(super) second: u32,
    #[pyo3(get)]
    pub(super) min_distance: f32,
}

#[pyclass(name = "ContactMap", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyContactMap {
    #[pyo3(get)]
    pub(super) residue_count: usize,
    #[pyo3(get)]
    pub(super) contacts: Vec<PyResidueContact>,
}

#[pyclass(name = "HalfSphereExposure", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyHalfSphereExposure {
    #[pyo3(get)]
    pub(super) residue: u32,
    #[pyo3(get)]
    pub(super) upper: u32,
    #[pyo3(get)]
    pub(super) lower: u32,
}

#[pyclass(name = "NativeContacts", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyNativeContacts {
    #[pyo3(get)]
    pub(super) native: usize,
    #[pyo3(get)]
    pub(super) kept: usize,
    #[pyo3(get)]
    pub(super) fraction: f64,
}

#[pyclass(name = "SurfaceContactOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySurfaceContactOptions {
    tolerance: f32,
    probe: f32,
    surface_density: f32,
    minimum_area: f32,
    backend: PySpatialBackend,
}

#[pymethods]
impl PySurfaceContactOptions {
    #[new]
    fn new(
        tolerance: f32,
        probe: f32,
        surface_density: f32,
        minimum_area: f32,
        backend: PySpatialBackend,
    ) -> Self {
        Self {
            tolerance,
            probe,
            surface_density,
            minimum_area,
            backend,
        }
    }
}

pub(crate) fn surface_contacts_analysis(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::analysis::Contact>>,
) -> PyResult<PyAnalysis> {
    super::contact_analysis_to_py(py, analysis)
}

impl PySurfaceContactOptions {
    pub(crate) fn native_parts(&self) -> (f32, f32, f32, f32, pdbiox::SpatialBackend) {
        (
            self.tolerance,
            self.probe,
            self.surface_density,
            self.minimum_area,
            self.backend.into(),
        )
    }

    pub(crate) fn from_native_parts(
        tolerance: f32,
        probe: f32,
        surface_density: f32,
        minimum_area: f32,
        backend: pdbiox::SpatialBackend,
    ) -> Self {
        Self {
            tolerance,
            probe,
            surface_density,
            minimum_area,
            backend: backend.into(),
        }
    }
}

#[pymethods]
impl PyStructure {
    #[pyo3(signature = (dictionary, options, context=None))]
    fn base_pairs(
        &self,
        py: Python<'_>,
        dictionary: &PyComponentDictionary,
        options: PyBasePairOptions,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<PyBasePair>> {
        let structure = self.structure().clone();
        let dictionary = dictionary.0.clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::base_pairs(&structure, dictionary.as_ref(), options.0, &context)
        })
        .map(|values| values.into_iter().map(PyBasePair::from).collect())
        .map_err(value_error)
    }

    #[pyo3(signature = (cutoff, minimum_separation, backend, context=None))]
    fn residue_contact_map(
        &self,
        py: Python<'_>,
        cutoff: f32,
        minimum_separation: u32,
        backend: PySpatialBackend,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<PyContactMap> {
        let structure = self.structure().clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::residue_contact_map(
                &structure,
                cutoff,
                minimum_separation,
                backend.into(),
                &context,
            )
        })
        .map(PyContactMap::from)
        .map_err(value_error)
    }

    #[pyo3(signature = (radius, backend, context=None))]
    fn half_sphere_exposure(
        &self,
        py: Python<'_>,
        radius: f32,
        backend: PySpatialBackend,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<PyHalfSphereExposure>> {
        let structure = self.structure().clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::half_sphere_exposure(&structure, radius, backend.into(), &context)
        })
        .map(|values| values.into_iter().map(PyHalfSphereExposure::from).collect())
        .map_err(value_error)
    }

    #[pyo3(signature = (first, second, cutoff, backend, context=None))]
    fn chain_interface(
        &self,
        py: Python<'_>,
        first: String,
        second: String,
        cutoff: f32,
        backend: PySpatialBackend,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<Vec<u32>> {
        let structure = self.structure().clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::chain_interface(
                &structure,
                &first,
                &second,
                cutoff,
                backend.into(),
                &context,
            )
        })
        .map(|values| values.into_iter().map(pdbiox::ResidueIndex::get).collect())
        .map_err(value_error)
    }

    #[pyo3(signature = (target, cutoff, tolerance, backend, context=None))]
    fn native_contact_fraction(
        &self,
        py: Python<'_>,
        target: &PyStructure,
        cutoff: f32,
        tolerance: f32,
        backend: PySpatialBackend,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<PyNativeContacts> {
        let reference = self.structure().clone();
        let target = target.structure().clone();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::native_contact_fraction(
                &reference,
                &target,
                cutoff,
                tolerance,
                backend.into(),
                &context,
            )
        })
        .map(PyNativeContacts::from)
        .map_err(value_error)
    }

    #[pyo3(signature = (radii, options, context=None))]
    fn surface_contacts(
        &self,
        py: Python<'_>,
        radii: PyReadonlyArray1<'_, f32>,
        options: PySurfaceContactOptions,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<super::PyContactTable> {
        let structure = self.structure().clone();
        let radius_values = radii.as_slice()?.to_vec();
        let context = context.map_or_else(default_context, PyExecutionContext::native);
        py.detach(move || {
            pdbiox::analysis::surface_contacts(
                &structure,
                &radius_values,
                pdbiox::analysis::SurfaceContactOptions {
                    tolerance: options.tolerance,
                    probe: options.probe,
                    surface_density: options.surface_density,
                    minimum_area: options.minimum_area,
                    backend: options.backend.into(),
                },
                &context,
            )
            .map(super::PyContactTable::from)
        })
        .map_err(value_error)
    }
}

#[pyfunction]
#[pyo3(signature = (structure, dictionary, options, context=None))]
pub(crate) fn base_pairs(
    py: Python<'_>,
    structure: &PyStructure,
    dictionary: &PyComponentDictionary,
    options: PyBasePairOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<PyBasePair>> {
    let structure = structure.structure().clone();
    let dictionary = dictionary.0.clone();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::base_pairs(&structure, dictionary.as_ref(), options.0, &context)
    })
    .map(|values| values.into_iter().map(PyBasePair::from).collect())
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, cutoff, minimum_separation, backend, context=None))]
pub(crate) fn residue_contact_map(
    py: Python<'_>,
    structure: &PyStructure,
    cutoff: f32,
    minimum_separation: u32,
    backend: PySpatialBackend,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyContactMap> {
    let structure = structure.structure().clone();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::residue_contact_map(
            &structure,
            cutoff,
            minimum_separation,
            backend.into(),
            &context,
        )
    })
    .map(PyContactMap::from)
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, radius, backend, context=None))]
pub(crate) fn half_sphere_exposure(
    py: Python<'_>,
    structure: &PyStructure,
    radius: f32,
    backend: PySpatialBackend,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<PyHalfSphereExposure>> {
    let structure = structure.structure().clone();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::half_sphere_exposure(&structure, radius, backend.into(), &context)
    })
    .map(|values| values.into_iter().map(PyHalfSphereExposure::from).collect())
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, first, second, cutoff, backend, context=None))]
pub(crate) fn chain_interface(
    py: Python<'_>,
    structure: &PyStructure,
    first: &str,
    second: &str,
    cutoff: f32,
    backend: PySpatialBackend,
    context: Option<&PyExecutionContext>,
) -> PyResult<Vec<u32>> {
    let structure = structure.structure().clone();
    let first = first.to_owned();
    let second = second.to_owned();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::chain_interface(
            &structure,
            &first,
            &second,
            cutoff,
            backend.into(),
            &context,
        )
    })
    .map(|values| values.into_iter().map(pdbiox::ResidueIndex::get).collect())
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (reference, target, cutoff, tolerance, backend, context=None))]
pub(crate) fn native_contact_fraction(
    py: Python<'_>,
    reference: &PyStructure,
    target: &PyStructure,
    cutoff: f32,
    tolerance: f32,
    backend: PySpatialBackend,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyNativeContacts> {
    let reference = reference.structure().clone();
    let target = target.structure().clone();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::native_contact_fraction(
            &reference,
            &target,
            cutoff,
            tolerance,
            backend.into(),
            &context,
        )
    })
    .map(PyNativeContacts::from)
    .map_err(value_error)
}

#[pyfunction]
#[pyo3(signature = (structure, radii, options, context=None))]
pub(crate) fn surface_contacts(
    py: Python<'_>,
    structure: &PyStructure,
    radii: PyReadonlyArray1<'_, f32>,
    options: PySurfaceContactOptions,
    context: Option<&PyExecutionContext>,
) -> PyResult<super::PyContactTable> {
    let structure = structure.structure().clone();
    let radii = radii.as_slice()?.to_vec();
    let context = context.map_or_else(default_context, PyExecutionContext::native);
    py.detach(move || {
        pdbiox::analysis::surface_contacts(
            &structure,
            &radii,
            pdbiox::analysis::SurfaceContactOptions {
                tolerance: options.tolerance,
                probe: options.probe,
                surface_density: options.surface_density,
                minimum_area: options.minimum_area,
                backend: options.backend.into(),
            },
            &context,
        )
        .map(super::PyContactTable::from)
    })
    .map_err(value_error)
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
