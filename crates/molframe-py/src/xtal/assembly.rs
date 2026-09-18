//! Owned Python projections for biological assemblies.

use super::types::{PyAssemblyDef, PyOperator};
use crate::geometry::PyRigid;
use crate::graph::PySpatialBackend;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "AssemblySet", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAssemblySet(pub(crate) molframe::xtal::AssemblySet);

#[pymethods]
impl PyAssemblySet {
    fn __len__(&self) -> usize {
        self.0.len()
    }
    #[getter]
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    fn get(&self, id: &str) -> Option<PyAssemblyDef> {
        self.0.get(id).cloned().map(Into::into)
    }
    fn operator(&self, id: &str) -> Option<PyOperator> {
        self.0.operator(id).cloned().map(Into::into)
    }
    fn assemblies(&self) -> Vec<PyAssemblyDef> {
        self.0.assemblies().cloned().map(Into::into).collect()
    }
    fn operators(&self) -> Vec<PyOperator> {
        self.0.operators().cloned().map(Into::into).collect()
    }
}

#[pyclass(name = "ChainInstance", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyChainInstance {
    #[pyo3(get)]
    source_chain: u32,
    #[pyo3(get)]
    transform: PyRigid,
    #[pyo3(get)]
    instance_id: u32,
}

#[pyclass(name = "AtomInstance", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAtomInstance {
    #[pyo3(get)]
    source_atom: u32,
    #[pyo3(get)]
    transform: PyRigid,
    #[pyo3(get)]
    instance_id: u32,
}

#[pymethods]
impl PyAtomInstance {
    fn apply(&self, position: [f32; 3]) -> [f32; 3] {
        self.transform.0.apply(position)
    }
}

#[pyclass(name = "AssemblyNeighbor", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyAssemblyNeighbor {
    #[pyo3(get)]
    first_instance: u32,
    #[pyo3(get)]
    first_atom: u32,
    #[pyo3(get)]
    second_instance: u32,
    #[pyo3(get)]
    second_atom: u32,
    #[pyo3(get)]
    distance_squared: f32,
}

#[pyclass(name = "AssemblyView", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAssemblyView(pub(crate) molframe::xtal::AssemblyView);

#[pymethods]
impl PyAssemblyView {
    #[getter]
    fn id(&self) -> &str {
        self.0.id()
    }
    #[getter]
    fn instance_count(&self) -> usize {
        self.0.instance_count()
    }
    fn chains(&self) -> Vec<PyChainInstance> {
        self.0.chains().map(Into::into).collect()
    }
    fn atoms(&self) -> Vec<PyAtomInstance> {
        self.0.atoms().map(Into::into).collect()
    }
    #[pyo3(signature = (model=0))]
    fn positions(&self, model: usize) -> PyResult<Vec<(PyAtomInstance, Option<[f32; 3]>)>> {
        let model = model_index(model)?;
        Ok(self
            .0
            .positions(model)
            .map(|(instance, position)| (instance.into(), position))
            .collect())
    }
    fn neighbors(
        &self,
        model: usize,
        cutoff: f32,
        backend: PySpatialBackend,
    ) -> PyResult<Vec<PyAssemblyNeighbor>> {
        self.0
            .neighbors(
                model_index(model)?,
                cutoff,
                backend.into(),
                &crate::core::execution::default_context(),
            )
            .map(|values| values.into_iter().map(Into::into).collect())
            .map_err(value_error)
    }
}

#[pymethods]
impl PyStructure {
    fn assembly_set(&self) -> Option<PyAssemblySet> {
        molframe::xtal::AssemblyExt::assembly_set(self.structure())
            .cloned()
            .map(PyAssemblySet)
    }

    #[pyo3(signature = (id, *, limit=None))]
    fn assembly(&self, id: &str, limit: Option<usize>) -> PyResult<PyAssemblyView> {
        let result = match limit {
            Some(limit) => {
                let assemblies = molframe::xtal::AssemblyExt::assembly_set(self.structure())
                    .ok_or_else(|| value_error("structure has no assembly metadata"))?;
                molframe::xtal::AssemblyView::with_limit(self.structure(), assemblies, id, limit)
            }
            None => molframe::xtal::AssemblyExt::assembly(self.structure(), id),
        };
        result.map(PyAssemblyView).map_err(value_error)
    }
}

fn model_index(model: usize) -> PyResult<molframe::ModelIndex> {
    u32::try_from(model)
        .map(molframe::ModelIndex::new)
        .map_err(|_| value_error("model index exceeds u32"))
}

fn value_error(error: impl std::fmt::Display) -> pyo3::PyErr {
    pyo3::exceptions::PyValueError::new_err(error.to_string())
}

impl From<molframe::xtal::ChainInstance> for PyChainInstance {
    fn from(value: molframe::xtal::ChainInstance) -> Self {
        Self {
            source_chain: value.source_chain.get(),
            transform: PyRigid(value.transform),
            instance_id: value.instance_id.get(),
        }
    }
}

impl From<molframe::xtal::AtomInstance> for PyAtomInstance {
    fn from(value: molframe::xtal::AtomInstance) -> Self {
        Self {
            source_atom: value.source_atom.get(),
            transform: PyRigid(value.transform),
            instance_id: value.instance_id.get(),
        }
    }
}

impl From<molframe::xtal::AssemblyNeighbor> for PyAssemblyNeighbor {
    fn from(value: molframe::xtal::AssemblyNeighbor) -> Self {
        Self {
            first_instance: value.first_instance.get(),
            first_atom: value.first_atom.get(),
            second_instance: value.second_instance.get(),
            second_atom: value.second_atom.get(),
            distance_squared: value.distance_squared,
        }
    }
}
