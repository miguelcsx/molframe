//! Deterministic clustering and representatives over ensemble distance matrices.

mod agglomerative;
mod dbscan;
mod matrix;
mod medoid;
mod types;

pub use agglomerative::agglomerative_clustering;
pub use dbscan::dbscan_clustering;
pub use medoid::medoid;
pub use types::{Clustering, Linkage};

#[cfg(test)]
#[path = "clustering_tests.rs"]
mod tests;
