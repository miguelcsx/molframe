//! Transactional structure edits with native validation and ownership.

use crate::core_annotations::PyAtomAnnotation;
use crate::errors::read_error;
use crate::geometry::PyRigid;
use crate::query::PySelection;
use crate::structure::PyStructure;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

#[pyclass(name = "StructureEditor", skip_from_py_object)]
pub(crate) struct PyStructureEditor {
    inner: Option<pdbiox::StructureEditor>,
}

fn missing_editor() -> PyErr {
    PyRuntimeError::new_err("the structure editor has already been committed")
}

#[pymethods]
impl PyStructureEditor {
    fn set_annotation(
        &mut self,
        py: Python<'_>,
        name: &str,
        annotation: PyAtomAnnotation,
    ) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(missing_editor)?
            .set_annotation(name, annotation.into())
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }

    fn remove_annotation(&mut self, name: &str) -> PyResult<bool> {
        self.inner
            .as_mut()
            .ok_or_else(missing_editor)
            .map(|editor| editor.remove_annotation(name))
    }

    fn clear_extensions(&mut self) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(missing_editor)
            .map(pdbiox::StructureEditor::clear_extensions)
    }

    fn rename_chain(
        &mut self,
        py: Python<'_>,
        chain: &crate::index::PyChainIndex,
        label: &str,
    ) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(missing_editor)?
            .rename_chain(chain.0, label)
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }

    fn transform(
        &mut self,
        py: Python<'_>,
        selection: &PySelection,
        transform: &PyRigid,
    ) -> PyResult<()> {
        let selection = selection.inner.clone();
        let rigid = transform.0;
        self.inner
            .as_mut()
            .ok_or_else(missing_editor)?
            .transform(&selection, |position| rigid.apply(position))
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }

    fn delete_atoms(&mut self, py: Python<'_>, selection: &PySelection) -> PyResult<()> {
        self.inner
            .as_mut()
            .ok_or_else(missing_editor)?
            .delete_atoms(&selection.inner)
            .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
    }

    fn commit(&mut self, py: Python<'_>) -> PyResult<PyStructure> {
        let editor = self.inner.take().ok_or_else(missing_editor)?;
        py.detach(move || editor.commit())
            .map(PyStructure::new)
            .map_err(|findings| read_error(py, &findings))
    }
}

#[pymethods]
impl PyStructure {
    fn edit(&self) -> PyStructureEditor {
        PyStructureEditor {
            inner: Some(self.structure().edit()),
        }
    }

    fn materialize(&self, py: Python<'_>, selection: &PySelection) -> PyResult<PyStructure> {
        let structure = self.structure().clone();
        let selection = selection.inner.clone();
        py.detach(move || structure.materialize(&selection))
            .map(PyStructure::new)
            .map_err(|findings| read_error(py, &findings))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStructureEditor>()?;
    module.add("CoordinateEditor", module.getattr("CoordinateEdit")?)?;
    Ok(())
}
