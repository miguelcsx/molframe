//! Python conversion for values emitted by native declarative plans.

use super::{geometry, physical, spatial, structure, trajectory};
use crate::analysis::contact_analysis_to_py as contact_table_analysis_to_py;
use crate::contract::PyAnalysis;
use crate::structure::PyStructure;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyFloat};

pub(super) fn plan_result_dict<'py>(
    py: Python<'py>,
    native: pdbiox::PlanResult,
    structure: Option<&'py PyStructure>,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    let structure = structure.map(PyStructure::structure);
    for entry in native.entries {
        result.set_item(entry.id.as_ref(), plan_value(py, entry.value, structure)?)?;
    }
    Ok(result)
}

pub(super) fn contact_analysis_to_py(
    py: Python<'_>,
    analysis: pdbiox::Analysis<Vec<pdbiox::analysis::Contact>>,
) -> PyResult<PyAnalysis> {
    contact_table_analysis_to_py(py, analysis)
}

fn plan_value(
    py: Python<'_>,
    value: pdbiox::PlanValue,
    structure: Option<&pdbiox::Structure>,
) -> PyResult<Py<PyAny>> {
    match value {
        pdbiox::PlanValue::Contacts(analysis) => {
            Ok(Py::new(py, contact_analysis_to_py(py, *analysis)?)?.into_any())
        }
        pdbiox::PlanValue::Rmsd(value) => Ok(PyFloat::new(py, value).unbind().into_any()),
        pdbiox::PlanValue::Comparison(value) => {
            Ok(PyFloat::new(py, value.value).unbind().into_any())
        }
        pdbiox::PlanValue::Selection(value) => {
            Ok(Py::new(py, crate::query::PyEvaluation::from(*value))?.into_any())
        }
        pdbiox::PlanValue::Structure(value) => {
            Ok(Py::new(py, structure::value_to_python(py, *value, structure)?)?.into_any())
        }
        pdbiox::PlanValue::BondInference(report) => Ok(Py::new(
            py,
            crate::facade::PyBondInferenceReport::from_native(*report),
        )?
        .into_any()),
        pdbiox::PlanValue::Geometry(value) => geometry::value_to_python(py, *value),
        pdbiox::PlanValue::Spatial(value) => match *value {
            pdbiox::SpatialValue::NeighborPairs(values) => {
                Ok(Py::new(py, spatial::neighbor_table_from_pairs(py, values)?)?.into_any())
            }
            pdbiox::SpatialValue::AtomsWithin(value) => {
                Ok(spatial::atom_indices_to_python(py, value)
                    .unbind()
                    .into_any())
            }
        },
        pdbiox::PlanValue::Physical(value) => physical::value_to_python_public(py, *value)
            .map(|value| Py::new(py, value).map(pyo3::Py::into_any))?,
        pdbiox::PlanValue::Surface(value) => match *value {
            pdbiox::SurfaceValue::SolventAccessibleSurface(values) => {
                use numpy::IntoPyArray;
                Ok(numpy::ndarray::Array1::from_vec(values)
                    .into_pyarray(py)
                    .unbind()
                    .into_any())
            }
            pdbiox::SurfaceValue::BuriedSurface(value) => {
                Ok(Py::new(py, crate::surface_types::PyBuriedSurface::from(value))?.into_any())
            }
        },
        pdbiox::PlanValue::Trajectory(value) => trajectory::value_to_python(py, *value),
    }
}
