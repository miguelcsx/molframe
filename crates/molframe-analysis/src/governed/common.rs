use super::{AnalysisDescriptor, ForbiddenResolution, FrameKernelResult, Requirement};
use molframe_core::contract::{HydrogenPolicy, ParameterValue, PolicyField};
use molframe_core::structure::Structure;
use molframe_spatial::SpatialBackend;

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

/// An analysis whose geometry is the position of a modelled hydrogen: excluding the
/// hydrogens leaves it nothing to measure, and an input without them has no answer.
pub(super) fn needing_hydrogens(descriptor: AnalysisDescriptor) -> AnalysisDescriptor {
    descriptor
        .forbidding(ForbiddenResolution::new(
            PolicyField::Hydrogens,
            "the donor-hydrogen-acceptor geometry is the modelled hydrogen, and excluding \
             the hydrogens would return no bonds rather than an answer",
            |policy| matches!(policy.hydrogens, HydrogenPolicy::Exclude),
        ))
        .requiring(Requirement::ExplicitHydrogens)
}
