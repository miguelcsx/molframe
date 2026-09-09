//! Parsing and writing bindings for text trajectory records.
//!
//! Split from the record types so each file stays within the size ceiling.

use super::{
    AimsError, CharmmError, GroError, PyAimsAtom, PyAimsGeometry, PyCharmmAtom, PyCharmmCard,
    PyCharmmCardFormat, PyGroAtom, PyGroFrame, PyTxyzAtom, PyTxyzFrame, TxyzError,
};
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn parse_aims_geometry(py: Python<'_>, text: &str) -> PyResult<PyAimsGeometry> {
    py.detach(move || -> PyResult<PyAimsGeometry> {
        pdbiox::traj::parse_aims_geometry(text)
            .map(Into::into)
            .map_err(|error| AimsError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_aims_geometry(py: Python<'_>, geometry: PyAimsGeometry) -> String {
    py.detach(move || -> String { pdbiox::traj::write_aims_geometry(&geometry.into()) })
}

#[pyfunction]
pub(crate) fn parse_gro_records(py: Python<'_>, text: &str) -> PyResult<Vec<PyGroFrame>> {
    py.detach(move || -> PyResult<Vec<PyGroFrame>> {
        pdbiox::traj::parse_gro_records(text)
            .map(|frames| frames.into_iter().map(Into::into).collect())
            .map_err(|error| GroError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_gro(py: Python<'_>, frames: Vec<PyGroFrame>) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        let frames = frames.into_iter().map(Into::into).collect::<Vec<_>>();
        pdbiox::traj::write_gro(&frames).map_err(|error| GroError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn parse_txyz_records(py: Python<'_>, text: &str) -> PyResult<Vec<PyTxyzFrame>> {
    py.detach(move || -> PyResult<Vec<PyTxyzFrame>> {
        pdbiox::traj::parse_txyz_records(text)
            .map(|frames| frames.into_iter().map(Into::into).collect())
            .map_err(|error| TxyzError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_txyz(py: Python<'_>, frames: Vec<PyTxyzFrame>) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        let frames = frames.into_iter().map(Into::into).collect::<Vec<_>>();
        pdbiox::traj::write_txyz(&frames).map_err(|error| TxyzError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn parse_charmm_record(py: Python<'_>, text: &str) -> PyResult<PyCharmmCard> {
    py.detach(move || -> PyResult<PyCharmmCard> {
        pdbiox::traj::parse_charmm_record(text)
            .map(Into::into)
            .map_err(|error| CharmmError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_charmm_card(py: Python<'_>, card: PyCharmmCard) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        pdbiox::traj::write_charmm_card(&card.into())
            .map_err(|error| CharmmError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("AimsError", module.py().get_type::<AimsError>())?;
    module.add("GroError", module.py().get_type::<GroError>())?;
    module.add("TxyzError", module.py().get_type::<TxyzError>())?;
    module.add("CharmmError", module.py().get_type::<CharmmError>())?;
    module.add_class::<PyAimsAtom>()?;
    module.add_class::<PyAimsGeometry>()?;
    module.add_class::<PyGroAtom>()?;
    module.add_class::<PyGroFrame>()?;
    module.add_class::<PyTxyzAtom>()?;
    module.add_class::<PyTxyzFrame>()?;
    module.add_class::<PyCharmmCardFormat>()?;
    module.add_class::<PyCharmmAtom>()?;
    module.add_class::<PyCharmmCard>()?;
    module.add_function(wrap_pyfunction!(parse_aims_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(write_aims_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(parse_gro_records, module)?)?;
    module.add_function(wrap_pyfunction!(write_gro, module)?)?;
    module.add_function(wrap_pyfunction!(parse_txyz_records, module)?)?;
    module.add_function(wrap_pyfunction!(write_txyz, module)?)?;
    module.add_function(wrap_pyfunction!(parse_charmm_record, module)?)?;
    module.add_function(wrap_pyfunction!(write_charmm_card, module)?)?;
    Ok(())
}
