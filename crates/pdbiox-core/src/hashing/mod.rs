//! Hashing for keys the engine assigned itself.

mod mix;

pub use mix::{IdentityBuildHasher, IdentityHashMap, IdentityHashSet, IdentityHasher};
