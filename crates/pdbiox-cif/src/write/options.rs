//! Explicit decisions and failures for canonical CIF projection.

use std::fmt::{Display, Formatter};

/// Decisions that cannot be recovered from a lowered structure.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CifWriteOptions {
    block_id: Option<Box<str>>,
    generate_connection_ids: bool,
    connection_type_id: Option<Box<str>>,
}

impl CifWriteOptions {
    /// Starts with no invented identifiers.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            block_id: None,
            generate_connection_ids: false,
            connection_type_id: None,
        }
    }

    /// Uses this explicit data-block identifier instead of `_entry.id`.
    #[must_use]
    pub fn with_block_id(mut self, block_id: impl Into<Box<str>>) -> Self {
        self.block_id = Some(block_id.into());
        self
    }

    /// Allows deterministic `_struct_conn.id` values when connectivity exists.
    ///
    /// The core bond graph does not retain a deposited connection identifier,
    /// so canonical connectivity cannot be written unless this choice is made.
    #[must_use]
    pub const fn with_generated_connection_ids(mut self) -> Self {
        self.generate_connection_ids = true;
        self
    }

    /// Declares the `_struct_conn.conn_type_id` assigned to graph bonds.
    ///
    /// The generic bond graph retains order and provenance but not the source
    /// dictionary's connection classification, so it must be supplied.
    #[must_use]
    pub fn with_connection_type_id(mut self, value: impl Into<Box<str>>) -> Self {
        self.connection_type_id = Some(value.into());
        self
    }

    pub(crate) fn block_id(&self) -> Option<&str> {
        self.block_id.as_deref()
    }

    pub(crate) const fn generates_connection_ids(&self) -> bool {
        self.generate_connection_ids
    }

    pub(crate) fn connection_type_id(&self) -> Option<&str> {
        self.connection_type_id.as_deref()
    }
}

/// A missing or invalid value that prevents canonical output.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CifWriteError {
    /// Neither `_entry.id` nor an explicit block identifier is available.
    MissingBlockId,
    /// The explicit or stored data-block identifier is not a bare CIF token.
    InvalidBlockId(String),
    /// A model lacks its deposited model number.
    MissingModelNumber {
        /// Zero-based model position.
        model: u32,
    },
    /// A model position exceeds the core model-index representation.
    ModelIndexOverflow {
        /// Zero-based model position.
        model: usize,
    },
    /// A required atom-site item is unavailable.
    MissingAtomField {
        /// Zero-based atom position.
        atom: u32,
        /// Dictionary item name without the category prefix.
        field: &'static str,
    },
    /// An atom-site identifier is zero, which records an absent source value.
    MissingAtomSiteId {
        /// Zero-based atom position.
        atom: u32,
    },
    /// A bond points outside the atom table.
    InvalidBondAtom {
        /// Zero-based bond position.
        bond: usize,
        /// Referenced zero-based atom position.
        atom: u32,
    },
    /// A bond endpoint cannot be represented by label identifiers.
    MissingBondField {
        /// Zero-based bond position.
        bond: usize,
        /// One-based endpoint number.
        endpoint: u8,
        /// Dictionary item name without the category prefix.
        field: &'static str,
    },
    /// Connectivity needs identifiers that the lowered graph does not retain.
    ConnectionIdsNotEnabled,
    /// Connectivity has no explicitly declared dictionary connection type.
    MissingConnectionTypeId,
}

impl Display for CifWriteError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingBlockId => formatter.write_str(
                "canonical CIF needs `_entry.id` or an explicit `CifWriteOptions::with_block_id`",
            ),
            Self::InvalidBlockId(id) => write!(formatter, "`{id}` is not a valid CIF block id"),
            Self::MissingModelNumber { model } => {
                write!(
                    formatter,
                    "model at position {model} has no deposited number"
                )
            }
            Self::ModelIndexOverflow { model } => {
                write!(
                    formatter,
                    "model position {model} exceeds the supported index range"
                )
            }
            Self::MissingAtomField { atom, field } => {
                write!(formatter, "atom {atom} has no required `{field}`")
            }
            Self::MissingAtomSiteId { atom } => {
                write!(formatter, "atom {atom} has no source `_atom_site.id`")
            }
            Self::InvalidBondAtom { bond, atom } => {
                write!(formatter, "bond {bond} refers to absent atom {atom}")
            }
            Self::MissingBondField {
                bond,
                endpoint,
                field,
            } => write!(
                formatter,
                "bond {bond} endpoint {endpoint} has no required `{field}`"
            ),
            Self::ConnectionIdsNotEnabled => formatter.write_str(
                "canonical connectivity requires explicit generated connection identifiers",
            ),
            Self::MissingConnectionTypeId => formatter.write_str(
                "canonical connectivity requires an explicit `_struct_conn.conn_type_id`",
            ),
        }
    }
}

impl std::error::Error for CifWriteError {}

pub(crate) fn valid_block_id(id: &str) -> bool {
    !id.is_empty()
        && id.bytes().all(|byte| {
            byte.is_ascii_graphic() && !matches!(byte, b'#' | b'\'' | b'"' | b'[' | b']')
        })
}
