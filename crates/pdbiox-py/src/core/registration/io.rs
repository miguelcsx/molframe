//! Registration for the public format-dispatch and variant I/O functions.

use crate::core_io::write_with_options;
use crate::io::{
    read_bytes, read_mmtf, read_pdb, read_pdbqt, read_pqr, read_with_diagnostics,
    read_with_options, write, write_bcif, write_bcif_with_options, write_mmcif, write_mmtf,
    write_pdb, write_pdbqt, write_pqr,
};
use pyo3::prelude::*;

pub(super) fn register_io_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(read_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(read_with_diagnostics, module)?)?;
    module.add_function(wrap_pyfunction!(read_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(read_mmtf, module)?)?;
    module.add_function(wrap_pyfunction!(read_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(read_pqr, module)?)?;
    module.add_function(wrap_pyfunction!(read_pdbqt, module)?)?;
    module.add_function(wrap_pyfunction!(write, module)?)?;
    module.add_function(wrap_pyfunction!(write_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(write_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(write_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(write_bcif_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(write_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(write_mmtf, module)?)?;
    module.add_function(wrap_pyfunction!(write_pqr, module)?)?;
    module.add_function(wrap_pyfunction!(write_pdbqt, module)?)?;
    Ok(())
}
