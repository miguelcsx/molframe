//! Deterministic clustering and representatives over ensemble distance matrices.

use crate::numeric::f64_from_usize;
use crate::{EnsembleDistanceMatrix, EnsembleGeometryError};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

include!("clustering/types.rs");
include!("clustering/agglomerative.rs");
include!("clustering/dbscan.rs");
include!("clustering/medoid.rs");
include!("clustering/helpers.rs");

#[cfg(test)]
#[path = "clustering_tests.rs"]
mod tests;
