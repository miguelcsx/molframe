//! `PyO3` classes for geometry operation nodes.

use super::super::dispatch::require_operation_tag;
use super::model::PyGeometryOperation;
use super::{
    eigen_from_dict, execute_operation, explain_operation, operation_to_dict, position_from_dict,
    principal_from_dict, validate_operation, weighted_from_dict,
};
use crate::geometry::PyEigenOptions;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

macro_rules! position_operation {
    ($type:ident, $python:literal, $tag:literal, $variant:ident) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $type {
            pub(crate) operation: PyGeometryOperation,
        }

        #[pymethods]
        impl $type {
            #[new]
            #[pyo3(signature = (*, positions))]
            fn new(py: Python<'_>, positions: Py<PyAny>) -> PyResult<Self> {
                let operation = PyGeometryOperation::$variant { positions };
                validate_operation(py, &operation)?;
                Ok(Self { operation })
            }

            fn execute(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                execute_operation(py, &self.operation)
            }

            fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                explain_operation(py, &self.operation)
            }

            fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                operation_to_dict(py, &self.operation)
            }

            #[staticmethod]
            fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
                Ok(Self {
                    operation: position_from_dict(config, $tag, |positions| {
                        PyGeometryOperation::$variant { positions }
                    })?,
                })
            }

            fn __repr__(&self) -> String {
                format!("{}(positions=<array>)", $python)
            }
        }
    };
}

position_operation!(PyCentroid, "Centroid", "centroid", Centroid);
position_operation!(
    PyDistanceMatrixOperation,
    "DistanceMatrixOp",
    "distance_matrix",
    DistanceMatrix
);

macro_rules! weighted_operation {
    ($type:ident, $python:literal, $tag:literal, $variant:ident) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $type {
            pub(crate) operation: PyGeometryOperation,
        }

        #[pymethods]
        impl $type {
            #[new]
            #[pyo3(signature = (*, positions, masses=None))]
            fn new(
                py: Python<'_>,
                positions: Py<PyAny>,
                masses: Option<Py<PyAny>>,
            ) -> PyResult<Self> {
                let operation = PyGeometryOperation::$variant { positions, masses };
                validate_operation(py, &operation)?;
                Ok(Self { operation })
            }

            fn execute(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                execute_operation(py, &self.operation)
            }

            fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                explain_operation(py, &self.operation)
            }

            fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                operation_to_dict(py, &self.operation)
            }

            #[staticmethod]
            fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
                Ok(Self {
                    operation: weighted_from_dict(config, $tag, |positions, masses| {
                        PyGeometryOperation::$variant { positions, masses }
                    })?,
                })
            }

            fn __repr__(&self) -> String {
                format!("{}(positions=<array>, masses=<array-or-none>)", $python)
            }
        }
    };
}

weighted_operation!(
    PyCentreOfMass,
    "CentreOfMass",
    "centre_of_mass",
    CentreOfMass
);
weighted_operation!(
    PyRadiusOfGyration,
    "RadiusOfGyration",
    "radius_of_gyration",
    RadiusOfGyration
);
weighted_operation!(
    PyInertiaTensor,
    "InertiaTensor",
    "inertia_tensor",
    InertiaTensor
);

macro_rules! eigen_operation {
    ($type:ident, $python:literal, $tag:literal, $variant:ident) => {
        #[pyclass(name = $python, frozen)]
        pub(crate) struct $type {
            pub(crate) operation: PyGeometryOperation,
        }

        #[pymethods]
        impl $type {
            #[new]
            #[pyo3(signature = (*, positions, options=None))]
            fn new(
                py: Python<'_>,
                positions: Py<PyAny>,
                options: Option<PyEigenOptions>,
            ) -> PyResult<Self> {
                let operation = PyGeometryOperation::$variant {
                    positions,
                    options: options
                        .map_or_else(pdbiox::EigenOptions::standard, |value| value.inner),
                };
                validate_operation(py, &operation)?;
                Ok(Self { operation })
            }

            fn execute(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                execute_operation(py, &self.operation)
            }

            fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                explain_operation(py, &self.operation)
            }

            fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
                operation_to_dict(py, &self.operation)
            }

            #[staticmethod]
            fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
                Ok(Self {
                    operation: eigen_from_dict(config, $tag, |positions, options| {
                        PyGeometryOperation::$variant { positions, options }
                    })?,
                })
            }

            fn __repr__(&self) -> String {
                format!("{}(positions=<array>, options=<EigenOptions>)", $python)
            }
        }
    };
}

eigen_operation!(PyAsphericity, "Asphericity", "asphericity", Asphericity);
eigen_operation!(
    PyGyrationAxes,
    "GyrationAxes",
    "gyration_axes",
    GyrationAxes
);

#[pyclass(name = "PrincipalAxes", frozen)]
pub(crate) struct PyPrincipalAxes {
    pub(crate) operation: PyGeometryOperation,
}

#[pymethods]
impl PyPrincipalAxes {
    #[new]
    #[pyo3(signature = (*, positions, masses=None, options=None))]
    fn new(
        py: Python<'_>,
        positions: Py<PyAny>,
        masses: Option<Py<PyAny>>,
        options: Option<PyEigenOptions>,
    ) -> PyResult<Self> {
        let operation = PyGeometryOperation::PrincipalAxes {
            positions,
            masses,
            options: options.map_or_else(pdbiox::EigenOptions::standard, |value| value.inner),
        };
        validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        Ok(Self {
            operation: principal_from_dict(config)?,
        })
    }

    fn __repr__(&self) -> String {
        "PrincipalAxes(positions=<array>, masses=<array-or-none>, options=<EigenOptions>)"
            .to_owned()
    }
}

#[pyclass(name = "DistanceMatrixBetween", frozen)]
pub(crate) struct PyDistanceMatrixBetween {
    pub(crate) operation: PyGeometryOperation,
}

#[pymethods]
impl PyDistanceMatrixBetween {
    #[new]
    #[pyo3(signature = (*, left, right))]
    fn new(py: Python<'_>, left: Py<PyAny>, right: Py<PyAny>) -> PyResult<Self> {
        let operation = PyGeometryOperation::DistanceMatrixBetween { left, right };
        validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, "distance_matrix_between")?;
        let operation = PyGeometryOperation::DistanceMatrixBetween {
            left: super::super::required(config, "left")?.unbind(),
            right: super::super::required(config, "right")?.unbind(),
        };
        validate_operation(config.py(), &operation)?;
        Ok(Self { operation })
    }

    fn __repr__(&self) -> String {
        "DistanceMatrixBetween(left=<array>, right=<array>)".to_owned()
    }
}

#[pyclass(name = "Rmsf", frozen)]
pub(crate) struct PyRmsf {
    pub(crate) operation: PyGeometryOperation,
}

#[pymethods]
impl PyRmsf {
    #[new]
    #[pyo3(signature = (*, frames))]
    fn new(py: Python<'_>, frames: Py<PyAny>) -> PyResult<Self> {
        let operation = PyGeometryOperation::Rmsf { frames };
        validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, "rmsf")?;
        let operation = PyGeometryOperation::Rmsf {
            frames: super::super::required(config, "frames")?.unbind(),
        };
        validate_operation(config.py(), &operation)?;
        Ok(Self { operation })
    }

    fn __repr__(&self) -> String {
        "Rmsf(frames=<array>)".to_owned()
    }
}

pub(super) fn typed_classes(value: &Bound<'_, PyAny>) -> PyResult<Option<PyGeometryOperation>> {
    macro_rules! extract {
        ($type:ty) => {
            if let Ok(value) = value.extract::<PyRef<'_, $type>>() {
                return Ok(Some(value.operation.clone_ref(value.py())));
            }
        };
    }
    extract!(PyCentroid);
    extract!(PyDistanceMatrixOperation);
    extract!(PyCentreOfMass);
    extract!(PyRadiusOfGyration);
    extract!(PyInertiaTensor);
    extract!(PyPrincipalAxes);
    extract!(PyAsphericity);
    extract!(PyGyrationAxes);
    extract!(PyDistanceMatrixBetween);
    extract!(PyRmsf);
    Ok(None)
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCentroid>()?;
    module.add_class::<PyDistanceMatrixOperation>()?;
    module.add_class::<PyCentreOfMass>()?;
    module.add_class::<PyRadiusOfGyration>()?;
    module.add_class::<PyInertiaTensor>()?;
    module.add_class::<PyPrincipalAxes>()?;
    module.add_class::<PyAsphericity>()?;
    module.add_class::<PyGyrationAxes>()?;
    module.add_class::<PyDistanceMatrixBetween>()?;
    module.add_class::<PyRmsf>()?;
    Ok(())
}

#[cfg(test)]
#[path = "classes_tests.rs"]
mod tests;
