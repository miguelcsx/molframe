//! The one error the comparison scores share.

use pdbiox_geom::SuperposeError;

/// Why a comparison could not be scored.
#[derive(Debug, thiserror::Error)]
pub enum CompareError {
    /// The two coordinate sets differ in length, so their indices cannot name
    /// the same points.
    #[error("model has {model} points but reference has {reference}")]
    LengthMismatch {
        /// Number of points in the model.
        model: usize,
        /// Number of points in the reference.
        reference: usize,
    },
    /// The superposition a score depends on could not be found.
    #[error("superposition failed: {0:?}")]
    Superpose(SuperposeError),
    /// Exact chemical mapping exceeded its explicit enumeration bound.
    #[error("atom mapping exceeded the limit of {limit} alternatives")]
    MappingLimit {
        /// Caller-selected upper bound.
        limit: usize,
    },
    /// A component without atoms cannot define an RMSD.
    #[error("component has no atoms to compare")]
    EmptyComponent,
    /// A correspondence is empty, out of bounds or not one-to-one.
    #[error("comparison mapping is empty, out of bounds or not one-to-one")]
    InvalidMapping,
    /// A distance cutoff is not finite and strictly positive.
    #[error("comparison distance cutoff must be finite and greater than zero")]
    InvalidDistanceCutoff,
    /// A score's length scale is not finite and strictly positive.
    #[error("comparison length scales must be finite and greater than zero")]
    InvalidLengthScale,
    /// A comparison received non-finite coordinates or score parameters.
    #[error("comparison coordinates and score parameters must be finite")]
    InvalidScoreInput,
    /// No pair met the caller's comparison-domain criteria.
    #[error("comparison domain contains no eligible pairs")]
    NoComparablePairs,
    /// An unqualified chain name was supplied under an explicit namespace policy.
    #[error("comparison requires a concrete label or auth identifier namespace")]
    UnsupportedNamespace,
}
