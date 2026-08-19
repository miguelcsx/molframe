//! Typed, reusable facade operations.
//!
//! Requests and results are defined separately from the deterministic executor
//! so each capability family can add native variants without creating a
//! stringly-typed dispatcher or a monolithic facade module.

#[cfg(feature = "compare")]
mod comparison;
mod geometry;
mod physical;
mod plan;
mod requests;
mod spatial;
mod spatial_cache;
mod structure;
#[cfg(feature = "surface")]
mod surface;
#[cfg(feature = "traj")]
mod trajectory;

#[cfg(feature = "compare")]
pub use comparison::{ComparisonMetric, ComparisonRequest, ComparisonResult};
pub use geometry::{GeometryRequest, GeometryValue};
pub use physical::{PhysicalRequest, PhysicalValue};
pub use plan::Plan;
pub use requests::FloatInput;
#[cfg(feature = "surface")]
pub use requests::MaskInput;
pub use requests::{
    ContactsRequest, CoordinateInput, ExecutionPlanError, FrameInput, IndexInput, PlanInput,
    PlanOperation, PlanResult, PlanResultEntry, PlanValue, RmsdRequest, ScalarInput,
    SelectionRequest,
};
pub use spatial::{SpatialRequest, SpatialValue};
pub use structure::{StructureRequest, StructureValue};
#[cfg(feature = "surface")]
pub use surface::{SurfaceRequest, SurfaceValue};
#[cfg(feature = "traj")]
pub use trajectory::{TrajectoryRequest, TrajectoryValue};

#[cfg(test)]
#[path = "operations_tests.rs"]
mod tests;
