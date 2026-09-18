//! Immutable metadata configuration for reciprocal-space tables.

use super::reflection::PyReflectionDataset;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "ReflectionMetadata", frozen, from_py_object)]
#[derive(Clone, Debug, Default)]
pub(crate) struct PyReflectionMetadata {
    pub(crate) title: Box<str>,
    pub(crate) cell: Option<molframe::UnitCell>,
    pub(crate) space_group_number: Option<i32>,
    pub(crate) space_group_name: Option<Box<str>>,
    pub(crate) datasets: Vec<molframe::xtal::ReflectionDataset>,
    pub(crate) history: Vec<Box<str>>,
    pub(crate) symmetry_operations: Vec<Box<str>>,
    pub(crate) sort_order: [i32; 5],
    pub(crate) resolution_range: Option<[f64; 2]>,
    pub(crate) missing_value: Option<f32>,
    pub(crate) extra_header_records: Vec<Box<str>>,
}

#[pymethods]
impl PyReflectionMetadata {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    fn with_title(&self, title: String) -> Self {
        Self {
            title: title.into(),
            ..self.clone()
        }
    }

    fn with_crystal(
        &self,
        cell: ([f64; 3], [f64; 3]),
        space_group_number: i32,
        space_group_name: String,
    ) -> PyResult<Self> {
        if space_group_number <= 0 || space_group_name.trim().is_empty() {
            return Err(PyValueError::new_err(
                "space_group_number must be positive and space_group_name non-empty",
            ));
        }
        Ok(Self {
            cell: Some(molframe::UnitCell {
                lengths: cell.0,
                angles: cell.1,
            }),
            space_group_number: Some(space_group_number),
            space_group_name: Some(space_group_name.into()),
            ..self.clone()
        })
    }

    fn with_symmetry_operations(&self, symmetry_operations: Vec<String>) -> PyResult<Self> {
        if symmetry_operations.is_empty()
            || symmetry_operations
                .iter()
                .any(|operation| operation.trim().is_empty())
        {
            return Err(PyValueError::new_err(
                "symmetry_operations must contain non-empty triplets",
            ));
        }
        Ok(Self {
            symmetry_operations: symmetry_operations.into_iter().map(Into::into).collect(),
            ..self.clone()
        })
    }

    fn with_datasets(&self, datasets: &Bound<'_, PyAny>) -> PyResult<Self> {
        let mut native = Vec::new();
        for dataset in datasets.try_iter()? {
            native.push(
                dataset?
                    .extract::<PyRef<'_, PyReflectionDataset>>()?
                    .0
                    .clone(),
            );
        }
        Ok(Self {
            datasets: native,
            ..self.clone()
        })
    }

    fn with_history(&self, history: Vec<String>) -> Self {
        Self {
            history: history.into_iter().map(Into::into).collect(),
            ..self.clone()
        }
    }

    fn with_sort_order(&self, sort_order: [i32; 5]) -> Self {
        Self {
            sort_order,
            ..self.clone()
        }
    }

    fn with_resolution_range(&self, resolution_range: [f64; 2]) -> Self {
        Self {
            resolution_range: Some(resolution_range),
            ..self.clone()
        }
    }

    fn with_missing_value(&self, missing_value: f32) -> Self {
        Self {
            missing_value: Some(missing_value),
            ..self.clone()
        }
    }

    fn with_extra_header_records(&self, extra_header_records: Vec<String>) -> Self {
        Self {
            extra_header_records: extra_header_records.into_iter().map(Into::into).collect(),
            ..self.clone()
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ReflectionMetadata(title={:?}, has_crystal={}, symmetry_operations={})",
            self.title,
            self.cell.is_some(),
            self.symmetry_operations.len()
        )
    }
}
