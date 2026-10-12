//! Explicit decisions and failures for canonical CIF projection.

use std::fmt::{Display, Formatter};
use std::io;

/// Decisions that cannot be recovered from a lowered structure.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CifWriteOptions {
    block_id: Option<Box<str>>,
    generate_connection_ids: bool,
    extra_categories: Option<Box<str>>,
}

impl CifWriteOptions {
    /// Starts with deterministic connection identifiers and no explicit block id.
    ///
    /// A connection's identifier is not retained by the bond graph, so a
    /// canonical file numbers its connections in order. The numbering carries no
    /// meaning of its own; [`Self::without_generated_connection_ids`] refuses
    /// instead for a caller that wants nothing invented.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            block_id: None,
            generate_connection_ids: true,
            extra_categories: None,
        }
    }

    /// Refuses connectivity rather than numbering its connections.
    #[must_use]
    pub const fn without_generated_connection_ids(mut self) -> Self {
        self.generate_connection_ids = false;
        self
    }

    /// Appends already-rendered mmCIF categories to the output.
    ///
    /// For metadata the generic structure keeps in a typed extension (assemblies,
    /// non-crystallographic operators, space group) that this crate cannot see.
    /// The text must hold complete categories and no data-block header.
    #[must_use]
    pub fn with_extra_categories(mut self, text: impl Into<Box<str>>) -> Self {
        let text = text.into();
        self.extra_categories = (!text.is_empty()).then_some(text);
        self
    }

    /// Rendered categories appended by the caller.
    #[must_use]
    pub fn extra_categories(&self) -> Option<&str> {
        self.extra_categories.as_deref()
    }

    /// Uses this explicit data-block identifier instead of `_entry.id`.
    #[must_use]
    pub fn with_block_id(mut self, block_id: impl Into<Box<str>>) -> Self {
        self.block_id = Some(block_id.into());
        self
    }

    /// Allows deterministic `_struct_conn.id` values when connectivity exists.
    ///
    /// This is the default; it re-enables the choice after
    /// [`Self::without_generated_connection_ids`].
    /// The connection's own `conn_type_id` travels with the bond and is written
    /// from it, so no class is chosen here.
    #[must_use]
    pub const fn with_generated_connection_ids(mut self) -> Self {
        self.generate_connection_ids = true;
        self
    }

    /// The explicit data-block identifier, when one was chosen.
    #[must_use]
    pub fn block_id(&self) -> Option<&str> {
        self.block_id.as_deref()
    }

    pub(crate) const fn generates_connection_ids(&self) -> bool {
        self.generate_connection_ids
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
    /// The total number of atom-site rows exceeds the platform index space.
    AtomRowCountOverflow,
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
    /// A secondary range endpoint lacks a required label identifier.
    MissingSecondaryField {
        /// Zero-based residue position.
        residue: u32,
        /// Missing dictionary item.
        field: &'static str,
    },
    /// Secondary label identities are not unique in the projected topology.
    AmbiguousSecondaryIdentity {
        /// Zero-based residue position.
        residue: u32,
    },
    /// Global secondary categories cannot encode independent model annotations.
    ModelSpecificSecondary,
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
            Self::AtomRowCountOverflow => {
                formatter.write_str("canonical atom-site row count exceeds the platform range")
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
            Self::MissingSecondaryField { residue, field } => {
                write!(
                    formatter,
                    "secondary residue {residue} has no required `{field}`"
                )
            }
            Self::AmbiguousSecondaryIdentity { residue } => write!(
                formatter,
                "secondary residue {residue} has an ambiguous label identity"
            ),
            Self::ModelSpecificSecondary => formatter.write_str(
                "canonical secondary categories cannot encode independent model annotations",
            ),
            Self::ConnectionIdsNotEnabled => formatter.write_str(
                "canonical connectivity requires explicit generated connection identifiers",
            ),
        }
    }
}

impl std::error::Error for CifWriteError {}

/// A canonical projection or output-stream failure.
#[derive(Debug)]
pub enum CifWriteToError {
    /// The semantic structure cannot be represented canonically.
    Projection(CifWriteError),
    /// The destination stopped accepting bytes.
    Output(io::Error),
}

impl Display for CifWriteToError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Projection(error) => Display::fmt(error, formatter),
            Self::Output(error) => write!(formatter, "output stream failed: {error}"),
        }
    }
}

impl std::error::Error for CifWriteToError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Projection(error) => Some(error),
            Self::Output(error) => Some(error),
        }
    }
}

impl From<CifWriteError> for CifWriteToError {
    fn from(error: CifWriteError) -> Self {
        Self::Projection(error)
    }
}

pub(crate) fn valid_block_id(id: &str) -> bool {
    !id.is_empty()
        && id.bytes().all(|byte| {
            byte.is_ascii_graphic() && !matches!(byte, b'#' | b'\'' | b'"' | b'[' | b']')
        })
}
