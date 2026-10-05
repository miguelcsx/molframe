//! Editing a structure: staged topology changes and scoped coordinate rewrites.
//!
//! A structure is immutable, so an edit is a transaction that publishes a new
//! structure and leaves the original, and every view into it, untouched.

use crate::bindings::{PySelection, PyStructure};
use crate::execution::PyExecutionContext;
use molframe::{CoordinateEditor, ExecutionContext, ModelIndex};
use numpy::PyArray2;
use numpy::ndarray::ArrayView2;
use pyo3::prelude::*;
use std::sync::Mutex;

/// Staged topology edits, published together by `finish()`.
///
/// Deleting atoms and renaming chains are validated as a whole when the edit
/// is finished: either the new structure satisfies every structural invariant
/// or `finish()` raises and the original is unchanged.
#[derive(Debug)]
#[pyclass(name = "StructureEditor", skip_from_py_object, module = "molframe")]
pub(crate) struct PyStructureEditor {
    inner: Option<molframe::StructureEditor>,
}

impl PyStructureEditor {
    pub(crate) const fn new(editor: molframe::StructureEditor) -> Self {
        Self {
            inner: Some(editor),
        }
    }

    fn open(&mut self) -> PyResult<&mut molframe::StructureEditor> {
        self.inner
            .as_mut()
            .ok_or_else(|| crate::error::internal("editor has already been finished"))
    }
}

#[pymethods]
impl PyStructureEditor {
    /// Stages caller-authorised boolean roles; None explicitly means unavailable.
    fn set_boolean_annotation(&mut self, name: &str, values: Vec<Option<bool>>) -> PyResult<()> {
        use molframe::engine::core::Presence;
        use molframe::{AnnotationColumn, AtomAnnotation};
        let entries = values.into_iter().map(|value| match value {
            Some(value) => (value, Presence::Present),
            None => (false, Presence::Unknown),
        });
        let column = AnnotationColumn::from_entries(entries)
            .map_err(|error| crate::error::value(error.to_string()))?;
        self.open()?
            .set_annotation(name, AtomAnnotation::Boolean(column))
            .map_err(crate::error::kernel)
    }

    /// Stages a rename of the chain at hierarchy index `chain`.
    fn rename_chain(&mut self, chain: u32, label: &str) -> PyResult<()> {
        self.open()?
            .rename_chain(molframe::ChainIndex::new(chain), label)
            .map_err(crate::error::kernel)
    }

    /// Drops assemblies, symmetry and other extensions attached to the structure.
    ///
    /// Deleting atoms is refused (`MOLFRAME-E3014`) while they are attached,
    /// because they refer to atoms the edit may remove; clearing them is the
    /// explicit, lossy choice that allows it.
    fn clear_extensions(&mut self) -> PyResult<()> {
        self.open()?.clear_extensions();
        Ok(())
    }

    /// Stages the deletion of every atom `selection` covers.
    fn delete(&mut self, selection: &PySelection) -> PyResult<()> {
        self.open()?
            .delete(selection.native())
            .map_err(crate::error::kernel)
    }

    /// Starts a coordinate rewrite of the structure this edit was opened from.
    ///
    /// The copy of the coordinates is charged to `context`'s memory budget
    /// before anything is allocated.
    #[pyo3(signature = (*, context=None))]
    fn coordinates(
        &mut self,
        py: Python<'_>,
        context: Option<&PyExecutionContext>,
    ) -> PyResult<PyCoordinateEditor> {
        let editor = self.open()?;
        let governing = match context {
            Some(context) => context.build(context.token())?,
            None => ExecutionContext::default(),
        };
        let started = py
            .detach(|| editor.coordinates(&governing))
            .map_err(crate::error::kernel)?;
        Ok(PyCoordinateEditor {
            editor: Mutex::new(started),
            context: governing,
        })
    }

    /// Validates the staged changes and publishes the new structure.
    fn finish(&mut self) -> PyResult<PyStructure> {
        let Some(editor) = self.inner.take() else {
            return Err(crate::error::internal("editor has already been finished"));
        };
        editor.finish().map(PyStructure::new).map_err(|findings| {
            let message = findings
                .first()
                .map_or_else(String::new, |first| first.message().to_owned());
            crate::error::from_findings(&findings, &message)
        })
    }
}

/// A private, writable copy of a structure's coordinates.
///
/// `positions()` returns a `NumPy` view onto the copy; write to it, then
/// `commit()` publishes a new structure holding its own snapshot of the result.
/// The original structure is never touched, and a view stays valid for as long as
/// the editor does, so later writes cannot change a structure already published.
#[derive(Debug)]
#[pyclass(
    name = "CoordinateEditor",
    frozen,
    skip_from_py_object,
    module = "molframe"
)]
pub(crate) struct PyCoordinateEditor {
    editor: Mutex<CoordinateEditor>,
    /// Keeps the accounting the copy was charged to alive as long as the copy.
    context: ExecutionContext,
}

#[pymethods]
impl PyCoordinateEditor {
    /// The writable `(atoms, 3)` `float32` positions of one dense model.
    #[pyo3(signature = (model=0))]
    fn positions<'py>(slf: &Bound<'py, Self>, model: u32) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let (pointer, rows) = {
            let this = slf.get();
            let Ok(mut editor) = this.editor.lock() else {
                return Err(crate::error::internal("the editor was poisoned"));
            };
            let positions = editor
                .try_positions_mut(ModelIndex::new(model))
                .map_err(crate::error::kernel)?;
            (positions.as_mut_ptr().cast::<f32>(), positions.len())
        };
        // SAFETY: the buffer is owned by the editor, which `owner` keeps alive for
        // as long as the array exists; it is never reallocated, `[f32; 3]` is
        // contiguous, and nothing else in Rust reads it while Python holds the GIL.
        let view = unsafe { ArrayView2::from_shape_ptr((rows, 3), pointer) };
        // SAFETY: the NumPy base object is the editor that owns the allocation.
        let array = unsafe { PyArray2::borrow_from_array(&view, slf.clone().into_any()) };
        Ok(array)
    }

    /// Publishes a new structure with the coordinates as they are now.
    ///
    /// The structure receives its own copy, so editing on afterwards cannot
    /// change it; the editor stays usable.
    fn commit(&self) -> PyResult<PyStructure> {
        let Ok(editor) = self.editor.lock() else {
            return Err(crate::error::internal("the editor was poisoned"));
        };
        editor
            .snapshot_detached(&self.context)
            .map(|core| PyStructure::new(molframe::Structure::from(core)))
            .map_err(|findings| {
                let message = findings
                    .first()
                    .map_or_else(String::new, |first| first.message().to_owned());
                crate::error::from_findings(&findings, &message)
            })
    }
}
