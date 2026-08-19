use pyo3::prelude::*;

pub(crate) fn register_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    super::arrow_types::register_classes(module)?;
    super::dlpack_types::register_classes(module)
}
