//! Typed boundary adapters for external topology projections and downloads.

use crate::bonds::PyBondOrder;
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyType;

pyo3::create_exception!(pdbiox_adapters, DownloadError, PyRuntimeError);
pyo3::create_exception!(pdbiox_adapters, TopologyExportError, PyValueError);

#[pyclass(name = "DownloadOptions", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyDownloadOptions(pdbiox::adapters::DownloadOptions);

#[pymethods]
impl PyDownloadOptions {
    #[new]
    #[pyo3(signature = (*, max_bytes, timeout_seconds, redirect_limit=0))]
    fn new(max_bytes: u64, timeout_seconds: f64, redirect_limit: usize) -> PyResult<Self> {
        let timeout = std::time::Duration::try_from_secs_f64(timeout_seconds).map_err(|_| {
            PyValueError::new_err("timeout_seconds must be finite and non-negative")
        })?;
        Ok(Self(pdbiox::adapters::DownloadOptions {
            max_bytes,
            timeout,
            redirect_limit,
        }))
    }

    #[getter]
    fn max_bytes(&self) -> u64 {
        self.0.max_bytes
    }

    #[getter]
    fn timeout_seconds(&self) -> f64 {
        self.0.timeout.as_secs_f64()
    }

    #[getter]
    fn redirect_limit(&self) -> usize {
        self.0.redirect_limit
    }
}

#[pyclass(name = "VerifiedDownload", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyVerifiedDownload(pdbiox::adapters::VerifiedDownload);

#[pymethods]
impl PyVerifiedDownload {
    #[getter]
    fn bytes(&self) -> Vec<u8> {
        self.0.bytes.clone()
    }

    #[getter]
    fn sha256(&self) -> String {
        self.0.sha256.clone()
    }
}

#[pyfunction]
pub(crate) fn fetch_verified(
    py: Python<'_>,
    url: String,
    expected_sha256: String,
    options: &PyDownloadOptions,
) -> PyResult<PyVerifiedDownload> {
    py.detach(move || {
        pdbiox::adapters::fetch_verified(&url, &expected_sha256, options.0)
            .map(PyVerifiedDownload)
            .map_err(download_error)
    })
}

#[pyclass(name = "ExportChain", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyExportChain(pdbiox::adapters::ExportChain);

#[pymethods]
impl PyExportChain {
    #[getter]
    fn id(&self) -> String {
        self.0.id.clone()
    }

    #[getter]
    fn residues(&self) -> (usize, usize) {
        (self.0.residues.start, self.0.residues.end)
    }
}

#[pyclass(name = "ExportResidue", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyExportResidue(pdbiox::adapters::ExportResidue);

#[pymethods]
impl PyExportResidue {
    #[getter]
    fn name(&self) -> String {
        self.0.name.clone()
    }

    #[getter]
    fn number(&self) -> Option<i32> {
        self.0.number
    }

    #[getter]
    fn insertion_code(&self) -> Option<String> {
        self.0.insertion_code.clone()
    }

    #[getter]
    fn is_heterogen(&self) -> bool {
        self.0.is_heterogen
    }

    #[getter]
    fn chain(&self) -> usize {
        self.0.chain
    }

    #[getter]
    fn atoms(&self) -> (usize, usize) {
        (self.0.atoms.start, self.0.atoms.end)
    }
}

#[pyclass(name = "ExportAtom", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyExportAtom(pdbiox::adapters::ExportAtom);

#[pymethods]
impl PyExportAtom {
    #[getter]
    fn name(&self) -> String {
        self.0.name.clone()
    }

    #[getter]
    fn atomic_number(&self) -> u8 {
        self.0.atomic_number
    }

    #[getter]
    fn element_symbol(&self) -> String {
        self.0.element_symbol.clone()
    }

    #[getter]
    fn mass(&self) -> f64 {
        self.0.mass
    }

    #[getter]
    fn serial(&self) -> Option<u32> {
        self.0.serial
    }

    #[getter]
    fn formal_charge(&self) -> Option<i8> {
        self.0.formal_charge
    }

    #[getter]
    fn occupancy(&self) -> Option<f32> {
        self.0.occupancy
    }

    #[getter]
    fn b_factor(&self) -> Option<f32> {
        self.0.b_factor
    }

    #[getter]
    fn alternate_location(&self) -> Option<String> {
        self.0.alternate_location.clone()
    }

    #[getter]
    fn residue(&self) -> usize {
        self.0.residue
    }

    #[getter]
    fn position(&self) -> [f32; 3] {
        self.0.position
    }
}

#[pyclass(name = "ExportBond", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyExportBond(pdbiox::adapters::ExportBond);

#[pymethods]
impl PyExportBond {
    #[getter]
    fn atom_a(&self) -> usize {
        self.0.atom_a
    }

    #[getter]
    fn atom_b(&self) -> usize {
        self.0.atom_b
    }

    #[getter]
    fn order(&self) -> PyBondOrder {
        self.0.order.into()
    }
}

#[pyclass(name = "TopologyExport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTopologyExport(pdbiox::adapters::TopologyExport);

#[pymethods]
impl PyTopologyExport {
    #[classmethod]
    fn from_model(
        _class: &Bound<'_, PyType>,
        structure: &PyStructure,
        model: usize,
        namespace: PyNamespace,
    ) -> PyResult<Self> {
        let model = u32::try_from(model)
            .map(pdbiox::ModelIndex::new)
            .map_err(|_| PyValueError::new_err("model index exceeds the native index range"))?;
        pdbiox::adapters::TopologyExport::from_model(structure.structure(), model, namespace.into())
            .map(Self)
            .map_err(topology_error)
    }

    #[getter]
    fn chains(&self) -> Vec<PyExportChain> {
        self.0.chains.iter().cloned().map(PyExportChain).collect()
    }

    #[getter]
    fn residues(&self) -> Vec<PyExportResidue> {
        self.0
            .residues
            .iter()
            .cloned()
            .map(PyExportResidue)
            .collect()
    }

    #[getter]
    fn atoms(&self) -> Vec<PyExportAtom> {
        self.0.atoms.iter().cloned().map(PyExportAtom).collect()
    }

    #[getter]
    fn bonds(&self) -> Vec<PyExportBond> {
        self.0.bonds.iter().copied().map(PyExportBond).collect()
    }
}

fn download_error(error: pdbiox::adapters::DownloadError) -> PyErr {
    DownloadError::new_err(error.to_string())
}

fn topology_error(error: pdbiox::adapters::TopologyExportError) -> PyErr {
    TopologyExportError::new_err(error.to_string())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    module.add("DownloadError", py.get_type::<DownloadError>())?;
    module.add("TopologyExportError", py.get_type::<TopologyExportError>())?;
    module.add_class::<PyDownloadOptions>()?;
    module.add_class::<PyVerifiedDownload>()?;
    module.add_class::<PyExportChain>()?;
    module.add_class::<PyExportResidue>()?;
    module.add_class::<PyExportAtom>()?;
    module.add_class::<PyExportBond>()?;
    module.add_class::<PyTopologyExport>()?;
    module.add_function(wrap_pyfunction!(fetch_verified, module)?)?;
    Ok(())
}
