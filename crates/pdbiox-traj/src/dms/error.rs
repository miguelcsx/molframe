//! DMS failures with enough context to distinguish storage and model errors.

/// A storage, schema, or model error encountered while reading or writing DMS.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum DmsError {
    /// The destination could not be reserved without replacing another file.
    #[error("DMS filesystem operation failed: {0}")]
    Io(#[from] std::io::Error),
    /// `SQLite` could not open, read, or write the database.
    #[error("DMS SQLite operation failed: {0}")]
    Sqlite(#[from] rusqlite::Error),
    /// A required table is absent.
    #[error("DMS is missing required table {table}")]
    MissingTable {
        /// Required table name.
        table: &'static str,
    },
    /// A required column is absent.
    #[error("DMS table {table} is missing required column {column}")]
    MissingColumn {
        /// Table containing the required column.
        table: &'static str,
        /// Required column name.
        column: &'static str,
    },
    /// A standard column has an incompatible declared `SQLite` affinity.
    #[error("DMS column {table}.{column} has an incompatible declared type")]
    InvalidColumnType {
        /// Relation containing the column.
        table: &'static str,
        /// Standard DMS column name.
        column: &'static str,
    },
    /// Particle identifiers are not contiguous and zero-based.
    #[error("DMS particle ids must be contiguous and start at zero")]
    InvalidParticleIds,
    /// A particle record contains an invalid number or inconsistent vector.
    #[error("invalid DMS particle record with id {id}")]
    InvalidParticle {
        /// File particle identifier.
        id: i64,
    },
    /// A bond is reversed, duplicated, or references an absent particle.
    #[error("invalid DMS bond ({p0}, {p1})")]
    InvalidBond {
        /// First particle identifier.
        p0: i64,
        /// Second particle identifier.
        p1: i64,
    },
    /// A present global cell does not contain exactly vectors zero, one, and two.
    #[error("DMS global_cell must contain exactly ids 0, 1, and 2")]
    InvalidCell,
    /// An optional DMS version table is malformed.
    #[error("invalid DMS version table")]
    InvalidVersion,
    /// Topology and frame arrays have different particle counts.
    #[error("DMS topology and frame particle counts differ")]
    ParticleCountMismatch,
    /// The destination already exists; DMS writes never overwrite implicitly.
    #[error("DMS destination already exists")]
    DestinationExists,
}
