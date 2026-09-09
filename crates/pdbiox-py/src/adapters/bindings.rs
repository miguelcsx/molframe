//! Typed neutral boundaries for topology transfer and verified downloads.

use crate::bonds::PyBondOrder;
use crate::query::PyNamespace;
use crate::structure::PyStructure;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyType;

pyo3::create_exception!(pdbiox_adapters, DownloadError, PyRuntimeError);
pyo3::create_exception!(pdbiox_adapters, TopologyBatchError, PyValueError);

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

#[pyclass(name = "TopologyBatch", skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTopologyBatch(Option<pdbiox::adapters::TopologyBatch>);

#[pymethods]
impl PyTopologyBatch {
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
        pdbiox::adapters::TopologyBatch::from_model(structure.structure(), model, namespace.into())
            .map(|batch| Self(Some(batch)))
            .map_err(topology_error)
    }

    #[getter]
    fn strings(&self) -> PyResult<Vec<String>> {
        Ok(self.batch()?.strings.clone())
    }

    #[getter]
    fn chain_ids(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.chain_ids.clone())
    }

    #[getter]
    fn chain_residue_offsets(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.chain_residue_offsets.clone())
    }

    #[getter]
    fn residue_names(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.residue_names.clone())
    }

    #[getter]
    fn residue_numbers(&self) -> PyResult<Vec<i32>> {
        Ok(self.batch()?.residue_numbers.clone())
    }

    #[getter]
    fn residue_number_validity(&self) -> PyResult<Vec<bool>> {
        Ok(validity_flags(&self.batch()?.residue_number_validity))
    }

    #[getter]
    fn residue_insertion_codes(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.residue_insertion_codes.clone())
    }

    #[getter]
    fn residue_chain(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.residue_chain.clone())
    }

    #[getter]
    fn residue_is_heterogen(&self) -> PyResult<Vec<bool>> {
        let batch = self.batch()?;
        Ok((0..batch.residue_names.len())
            .map(|index| {
                u32::try_from(index).is_ok_and(|index| batch.residue_is_heterogen.test(index))
            })
            .collect())
    }

    #[getter]
    fn residue_atom_offsets(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.residue_atom_offsets.clone())
    }

    #[getter]
    fn atom_names(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.atom_names.clone())
    }

    #[getter]
    fn atomic_numbers(&self) -> PyResult<Vec<u8>> {
        Ok(self.batch()?.atomic_numbers.clone())
    }

    #[getter]
    fn masses(&self) -> PyResult<Vec<f64>> {
        Ok(self.batch()?.masses.clone())
    }

    #[getter]
    fn atom_serials(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.atom_serials.clone())
    }

    #[getter]
    fn atom_serial_validity(&self) -> PyResult<Vec<bool>> {
        Ok(validity_flags(&self.batch()?.atom_serial_validity))
    }

    #[getter]
    fn formal_charges(&self) -> PyResult<Vec<i8>> {
        Ok(self.batch()?.formal_charges.clone())
    }

    #[getter]
    fn formal_charge_validity(&self) -> PyResult<Vec<bool>> {
        Ok(validity_flags(&self.batch()?.formal_charge_validity))
    }

    #[getter]
    fn occupancies(&self) -> PyResult<Vec<f32>> {
        Ok(self.batch()?.occupancies.clone())
    }

    #[getter]
    fn occupancy_validity(&self) -> PyResult<Vec<bool>> {
        Ok(validity_flags(&self.batch()?.occupancy_validity))
    }

    #[getter]
    fn b_factors(&self) -> PyResult<Vec<f32>> {
        Ok(self.batch()?.b_factors.clone())
    }

    #[getter]
    fn b_factor_validity(&self) -> PyResult<Vec<bool>> {
        Ok(validity_flags(&self.batch()?.b_factor_validity))
    }

    #[getter]
    fn atom_alternate_locations(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.atom_alternate_locations.clone())
    }

    #[getter]
    fn atom_residue(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.atom_residue.clone())
    }

    #[getter]
    fn position_x(&self) -> PyResult<Vec<f32>> {
        Ok(self.batch()?.position_x.clone())
    }

    #[getter]
    fn position_y(&self) -> PyResult<Vec<f32>> {
        Ok(self.batch()?.position_y.clone())
    }

    #[getter]
    fn position_z(&self) -> PyResult<Vec<f32>> {
        Ok(self.batch()?.position_z.clone())
    }

    #[getter]
    fn bond_atom_a(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.bond_atom_a.clone())
    }

    #[getter]
    fn bond_atom_b(&self) -> PyResult<Vec<u32>> {
        Ok(self.batch()?.bond_atom_b.clone())
    }

    #[getter]
    fn bond_orders(&self) -> PyResult<Vec<PyBondOrder>> {
        Ok(self
            .batch()?
            .bond_orders
            .iter()
            .copied()
            .map(Into::into)
            .collect())
    }

    /// Transfers the projection into native columnar storage without cloning it.
    fn transfer_to_structure(&mut self) -> PyResult<PyStructure> {
        let batch = self.0.take().ok_or_else(consumed_batch)?;
        batch
            .into_structure()
            .map(PyStructure::new)
            .map_err(topology_import_error)
    }
}

impl PyTopologyBatch {
    fn batch(&self) -> PyResult<&pdbiox::adapters::TopologyBatch> {
        self.0.as_ref().ok_or_else(consumed_batch)
    }
}

fn validity_flags(validity: &pdbiox::core::ValidityMask) -> Vec<bool> {
    (0..validity.len())
        .map(|index| validity.get(index).is_present())
        .collect()
}

fn download_error(error: pdbiox::adapters::DownloadError) -> PyErr {
    DownloadError::new_err(error.to_string())
}

fn topology_error(error: pdbiox::adapters::TopologyBatchError) -> PyErr {
    TopologyBatchError::new_err(error.to_string())
}

fn topology_import_error(error: pdbiox::adapters::TopologyImportError) -> PyErr {
    TopologyBatchError::new_err(error.to_string())
}

fn consumed_batch() -> PyErr {
    PyRuntimeError::new_err("topology batch has already transferred its storage")
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    module.add("DownloadError", py.get_type::<DownloadError>())?;
    module.add("TopologyBatchError", py.get_type::<TopologyBatchError>())?;
    module.add_class::<PyDownloadOptions>()?;
    module.add_class::<PyVerifiedDownload>()?;
    module.add_class::<PyTopologyBatch>()?;
    module.add("MISSING_STRING", pdbiox::adapters::MISSING_STRING)?;
    module.add_function(wrap_pyfunction!(fetch_verified, module)?)?;
    Ok(())
}
