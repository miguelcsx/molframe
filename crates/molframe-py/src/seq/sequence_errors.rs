//! Stable Python exception classes for sequence facade errors.

use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, A2mError, PyValueError);
create_exception!(_native, AlignError, PyValueError);
create_exception!(_native, FastqError, PyValueError);
create_exception!(_native, MatrixError, PyValueError);
create_exception!(_native, MsaError, PyValueError);
create_exception!(_native, RegionAlignError, PyValueError);
create_exception!(_native, RegionError, PyValueError);
create_exception!(_native, SequenceFormatError, PyValueError);
create_exception!(_native, SimilarKmerError, PyValueError);

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("A2mError", module.py().get_type::<A2mError>())?;
    module.add("AlignError", module.py().get_type::<AlignError>())?;
    module.add("FastqError", module.py().get_type::<FastqError>())?;
    module.add("MatrixError", module.py().get_type::<MatrixError>())?;
    module.add("MsaError", module.py().get_type::<MsaError>())?;
    module.add(
        "RegionAlignError",
        module.py().get_type::<RegionAlignError>(),
    )?;
    module.add("RegionError", module.py().get_type::<RegionError>())?;
    module.add(
        "SequenceFormatError",
        module.py().get_type::<SequenceFormatError>(),
    )?;
    module.add(
        "SimilarKmerError",
        module.py().get_type::<SimilarKmerError>(),
    )?;
    Ok(())
}
