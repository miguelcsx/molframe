//! Registration of typed core-index values.

use crate::index::{
    PyAtomIndex, PyBondIndex, PyChainIndex, PyChunkId, PyEntityIndex, PyInstanceId, PyModelIndex,
    PyResidueIndex,
};
use pyo3::prelude::*;

pub(crate) fn register_indices(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAtomIndex>()?;
    module.add_class::<PyBondIndex>()?;
    module.add_class::<PyChainIndex>()?;
    module.add_class::<PyChunkId>()?;
    module.add_class::<PyEntityIndex>()?;
    module.add_class::<PyInstanceId>()?;
    module.add_class::<PyModelIndex>()?;
    module.add_class::<PyResidueIndex>()?;
    Ok(())
}
