//! Python sequence-index conversion at the language boundary.

use pyo3::prelude::*;

pub(crate) fn normalise_index(index: isize, count: usize) -> Option<usize> {
    let count = isize::try_from(count).ok()?;
    let position = if index < 0 {
        count.checked_add(index)?
    } else {
        index
    };
    if !(0..count).contains(&position) {
        return None;
    }
    usize::try_from(position).ok()
}

macro_rules! index_class {
    ($name:ident, $python:literal, $rust:path) => {
        #[pyclass(name = $python, frozen, eq, from_py_object)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) struct $name(pub(crate) $rust);

        #[pymethods]
        impl $name {
            #[new]
            fn new(value: u32) -> Self {
                Self(<$rust>::new(value))
            }

            #[staticmethod]
            fn from_raw(value: u32) -> Self {
                Self(<$rust>::new(value))
            }

            #[getter]
            fn value(&self) -> u32 {
                self.0.get()
            }

            fn get(&self) -> u32 {
                self.0.get()
            }

            fn as_usize(&self) -> usize {
                self.0.as_usize()
            }

            fn next(&self) -> Option<Self> {
                self.0.next().map(Self)
            }

            fn __int__(&self) -> u32 {
                self.0.get()
            }

            fn __index__(&self) -> u32 {
                self.0.get()
            }

            fn __repr__(&self) -> String {
                format!("{}({})", $python, self.0.get())
            }
        }
    };
}

index_class!(PyAtomIndex, "AtomIndex", molframe::AtomIndex);
index_class!(PyBondIndex, "BondIndex", molframe::BondIndex);
index_class!(PyChainIndex, "ChainIndex", molframe::ChainIndex);
index_class!(PyEntityIndex, "EntityIndex", molframe::EntityIndex);
index_class!(PyInstanceId, "InstanceId", molframe::InstanceId);
index_class!(PyModelIndex, "ModelIndex", molframe::ModelIndex);
index_class!(PyResidueIndex, "ResidueIndex", molframe::ResidueIndex);

macro_rules! global_id_class {
    ($name:ident, $python:literal, $rust:path) => {
        #[pyclass(name = $python, frozen, eq, from_py_object)]
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub(crate) struct $name(pub(crate) $rust);

        #[pymethods]
        impl $name {
            #[new]
            fn new(value: u64) -> Self {
                Self(<$rust>::new(value))
            }

            #[staticmethod]
            fn from_raw(value: u64) -> Self {
                Self(<$rust>::new(value))
            }

            #[getter]
            fn value(&self) -> u64 {
                self.0.get()
            }

            fn get(&self) -> u64 {
                self.0.get()
            }

            fn next(&self) -> Option<Self> {
                self.0.next().map(Self)
            }

            fn __int__(&self) -> u64 {
                self.0.get()
            }

            fn __index__(&self) -> u64 {
                self.0.get()
            }

            fn __repr__(&self) -> String {
                format!("{}({})", $python, self.0.get())
            }
        }
    };
}

global_id_class!(PyDatasetId, "DatasetId", molframe::DatasetId);
global_id_class!(PyChunkId, "ChunkId", molframe::ChunkId);
global_id_class!(PyLogicalRow, "LogicalRow", molframe::LogicalRow);

#[pyclass(name = "LocalRow", frozen, eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PyLocalRow(pub(crate) molframe::LocalRow);

#[pymethods]
impl PyLocalRow {
    #[new]
    fn new(value: u32) -> Self {
        Self(molframe::LocalRow::new(value))
    }

    #[getter]
    fn value(&self) -> u32 {
        self.0.get()
    }

    fn get(&self) -> u32 {
        self.0.get()
    }

    fn next(&self) -> Option<Self> {
        self.0.next().map(Self)
    }

    fn __int__(&self) -> u32 {
        self.0.get()
    }

    fn __index__(&self) -> u32 {
        self.0.get()
    }

    fn __repr__(&self) -> String {
        format!("LocalRow({})", self.0.get())
    }
}

#[cfg(test)]
#[path = "index_tests.rs"]
mod tests;
