//! Path-based trajectory reading.
//!
//! This module owns format selection and the shared bounded-input path. The two
//! siblings own the per-format decoders, split by whether the payload is
//! decoded as text or read as raw bytes.

mod binary;
mod text;

mod bytes;
mod materialized;
mod snapshot;
mod steps;

use bytes::read_bytes;
pub use materialized::read_trajectory_materialized;
use snapshot::read_snapshot;
use steps::{dcd_steps, frame_steps, integer_data};
