//! Python exception classes for functional-geometry stages.

use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, MotifError, PyValueError);
create_exception!(_native, MappingError, PyValueError);
create_exception!(_native, MeasurementError, PyValueError);
create_exception!(_native, EvaluationError, PyValueError);
create_exception!(_native, SpecificationError, PyValueError);
create_exception!(_native, AmeMeasurementError, PyValueError);

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("MotifError", module.py().get_type::<MotifError>())?;
    module.add("MappingError", module.py().get_type::<MappingError>())?;
    module.add(
        "MeasurementError",
        module.py().get_type::<MeasurementError>(),
    )?;
    module.add("EvaluationError", module.py().get_type::<EvaluationError>())?;
    module.add(
        "SpecificationError",
        module.py().get_type::<SpecificationError>(),
    )?;
    module.add(
        "AmeMeasurementError",
        module.py().get_type::<AmeMeasurementError>(),
    )?;
    Ok(())
}

pub(crate) fn motif_error(error: pdbiox::fx::MotifError) -> PyErr {
    MotifError::new_err(error.to_string())
}

pub(crate) fn mapping_error(error: pdbiox::fx::MappingError) -> PyErr {
    MappingError::new_err(error.to_string())
}

pub(crate) fn measurement_error(error: pdbiox::fx::MeasurementError) -> PyErr {
    MeasurementError::new_err(error.to_string())
}

pub(crate) fn evaluation_error(error: pdbiox::fx::EvaluationError) -> PyErr {
    EvaluationError::new_err(error.to_string())
}

pub(crate) fn specification_error(error: pdbiox::fx::SpecificationError) -> PyErr {
    SpecificationError::new_err(error.to_string())
}

pub(crate) fn ame_measurement_error(error: pdbiox::fx::AmeMeasurementError) -> PyErr {
    AmeMeasurementError::new_err(error.to_string())
}
