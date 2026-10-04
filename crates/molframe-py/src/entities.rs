//! The distinct species a structure contains, apart from the chains that copy them.

use crate::bindings::PyStructure;
use molframe::EntityIndex;
use pyo3::prelude::*;

/// One chemical species; several chains may be copies of it.
#[derive(Clone, Debug)]
#[pyclass(name = "Entity", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyEntity {
    parent: PyStructure,
    index: u32,
}

impl PyEntity {
    fn text(&self, symbol: Option<molframe::SymbolId>) -> Option<String> {
        symbol
            .and_then(|symbol| self.parent.inner.engine().resolve(symbol))
            .map(str::to_owned)
    }
}

#[pymethods]
impl PyEntity {
    #[getter]
    const fn index(&self) -> u32 {
        self.index
    }

    /// The identifier the file gave the entity.
    #[getter]
    fn id(&self) -> Option<String> {
        let table = &self.parent.inner.engine().data().topology.entities;
        self.text(table.id(EntityIndex::new(self.index)))
    }

    /// `"polymer"`, `"non_polymer"`, `"water"`, `"branched"` or `"unknown"`.
    #[getter]
    fn kind(&self) -> Option<String> {
        let table = &self.parent.inner.engine().data().topology.entities;
        table
            .kind(EntityIndex::new(self.index))
            .map(|kind| crate::policy::snake(kind.name()))
    }

    #[getter]
    fn description(&self) -> Option<String> {
        let table = &self.parent.inner.engine().data().topology.entities;
        self.text(table.description(EntityIndex::new(self.index)))
    }

    /// The sequence the entity should have, as component names. What a chain
    /// actually models may be shorter; the difference is the unmodelled region.
    #[getter]
    fn sequence(&self) -> Vec<String> {
        let table = &self.parent.inner.engine().data().topology.entities;
        table
            .canonical_sequence(EntityIndex::new(self.index))
            .iter()
            .map(|symbol| match self.text(Some(*symbol)) {
                Some(name) => name,
                None => String::new(),
            })
            .collect()
    }

    /// The chains that are copies of this entity.
    #[getter]
    fn chains(&self) -> Vec<u32> {
        self.parent
            .inner
            .chains()
            .iter()
            .filter(|chain| {
                chain
                    .entity()
                    .is_some_and(|entity| entity.get() == self.index)
            })
            .map(|chain| chain.index().get())
            .collect()
    }

    fn __repr__(&self) -> String {
        format!("Entity(index={}, kind={:?})", self.index, self.kind())
    }
}

/// The entities of a structure, in file order.
#[derive(Clone, Debug)]
#[pyclass(name = "Entities", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyEntities {
    parent: PyStructure,
}

impl PyEntities {
    pub(crate) const fn new(parent: PyStructure) -> Self {
        Self { parent }
    }

    fn count(&self) -> usize {
        self.parent.inner.engine().data().topology.entities.len()
    }
}

#[pymethods]
impl PyEntities {
    fn __len__(&self) -> usize {
        self.count()
    }

    fn __getitem__(&self, index: isize) -> PyResult<PyEntity> {
        let length = self.count();
        let resolved = if index < 0 {
            length.checked_sub(index.unsigned_abs())
        } else {
            Some(index.cast_unsigned())
        };
        match resolved {
            Some(position) if position < length => Ok(PyEntity {
                parent: self.parent.clone(),
                index: u32::try_from(position)
                    .map_err(|_| crate::error::index("entity index is out of range"))?,
            }),
            _ => Err(crate::error::index("entity index is out of range")),
        }
    }

    fn __repr__(&self) -> String {
        format!("Entities(count={})", self.count())
    }
}
