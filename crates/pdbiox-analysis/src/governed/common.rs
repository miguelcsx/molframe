use super::{AnalysisDescriptor, FrameKernelResult};
use pdbiox_core::contract::ParameterValue;
use pdbiox_core::structure::Structure;
use pdbiox_spatial::SpatialBackend;

pub(super) const ALGORITHM_VERSION: &str = "1";

pub(super) fn descriptor(name: &'static str) -> AnalysisDescriptor {
    AnalysisDescriptor::new(name, ALGORITHM_VERSION)
}

pub(super) fn float(value: impl Into<f64>) -> ParameterValue {
    ParameterValue::Float(value.into().to_bits())
}

pub(super) fn integer(value: usize) -> ParameterValue {
    match i64::try_from(value) {
        Ok(value) => ParameterValue::Integer(value),
        Err(_) => ParameterValue::Text(value.to_string().into()),
    }
}

pub(super) fn backend(value: SpatialBackend) -> ParameterValue {
    ParameterValue::Text(format!("{value:?}").into())
}

pub(super) fn complete<T>(structure: &Structure, value: T) -> FrameKernelResult<T> {
    FrameKernelResult::complete(value, structure.atom_count())
}
