//! Shared provenance parameters for governed comparison entry points.

use crate::CeOptions;
use molframe_core::contract::{Analysis, ParameterValue};

pub(crate) fn add_ce_parameters<T>(result: &mut Analysis<T>, options: CeOptions) {
    result.provenance = result
        .provenance
        .clone()
        .with_parameter("window_size", integer(options.window_size))
        .with_parameter("max_gap", integer(options.max_gap))
        .with_parameter("max_paths", integer(options.max_paths))
        .with_parameter("memory_limit_bytes", integer(options.memory_limit_bytes))
        .with_parameter(
            "fragment_similarity_threshold",
            ParameterValue::Float(options.fragment_similarity_threshold.to_bits()),
        )
        .with_parameter(
            "path_similarity_threshold",
            ParameterValue::Float(options.path_similarity_threshold.to_bits()),
        )
        .with_parameter(
            "significance",
            ParameterValue::Text(format!("{:?}", options.significance).into()),
        );
}

pub(crate) fn integer(value: usize) -> ParameterValue {
    match i64::try_from(value) {
        Ok(value) => ParameterValue::Integer(value),
        Err(_) => ParameterValue::Text(value.to_string().into()),
    }
}
