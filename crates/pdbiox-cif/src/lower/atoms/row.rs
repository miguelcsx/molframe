//! Borrowed coordinate rows shared by text and binary CIF readers.
//!
//! A reader exposes one row for the duration of [`AtomSiteRowSink::feed`]. The
//! lowering stage consumes every value synchronously and never retains the row
//! or any text borrowed from it. This keeps the seam compatible with both raw
//! text slices and dictionary-indexed binary columns without allocating a
//! document value for every cell.

use crate::parser::Rows;
use std::borrow::Cow;

/// One of the `atom_site` items the structure lowerer reads.
///
/// The set is closed: these are exactly the items a coordinate row can supply
/// to `Structure`, and a reader that carries others carries them for some other
/// consumer. Naming them as a type is what lets a row answer by index instead
/// of by comparing the caller's string against every item it holds — which, at
/// fifteen accessor calls per atom, was most of what lowering did.
#[doc(hidden)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum Field {
    /// `_atom_site.group_PDB`
    GroupPdb,
    /// `_atom_site.id`
    Id,
    /// `_atom_site.type_symbol`
    TypeSymbol,
    /// `_atom_site.label_atom_id`
    LabelAtomId,
    /// `_atom_site.auth_atom_id`
    AuthAtomId,
    /// `_atom_site.label_alt_id`
    LabelAltId,
    /// `_atom_site.label_comp_id`
    LabelCompId,
    /// `_atom_site.auth_comp_id`
    AuthCompId,
    /// `_atom_site.label_asym_id`
    LabelAsymId,
    /// `_atom_site.auth_asym_id`
    AuthAsymId,
    /// `_atom_site.label_entity_id`
    LabelEntityId,
    /// `_atom_site.label_seq_id`
    LabelSeqId,
    /// `_atom_site.auth_seq_id`
    AuthSeqId,
    /// `_atom_site.pdbx_PDB_ins_code`
    InsCode,
    /// `_atom_site.Cartn_x`
    CartnX,
    /// `_atom_site.Cartn_y`
    CartnY,
    /// `_atom_site.Cartn_z`
    CartnZ,
    /// `_atom_site.occupancy`
    Occupancy,
    /// `_atom_site.B_iso_or_equiv`
    BFactor,
    /// `_atom_site.pdbx_formal_charge`
    FormalCharge,
    /// `_atom_site.pdbx_PDB_model_num`
    ModelNum,
}

impl Field {
    /// How many fields there are, for a row that stores them by position.
    pub const COUNT: usize = 21;

    /// The item name this field carries in an `atom_site` category.
    #[must_use]
    pub const fn item(self) -> &'static str {
        match self {
            Self::GroupPdb => "group_PDB",
            Self::Id => "id",
            Self::TypeSymbol => "type_symbol",
            Self::LabelAtomId => "label_atom_id",
            Self::AuthAtomId => "auth_atom_id",
            Self::LabelAltId => "label_alt_id",
            Self::LabelCompId => "label_comp_id",
            Self::AuthCompId => "auth_comp_id",
            Self::LabelAsymId => "label_asym_id",
            Self::AuthAsymId => "auth_asym_id",
            Self::LabelEntityId => "label_entity_id",
            Self::LabelSeqId => "label_seq_id",
            Self::AuthSeqId => "auth_seq_id",
            Self::InsCode => "pdbx_PDB_ins_code",
            Self::CartnX => "Cartn_x",
            Self::CartnY => "Cartn_y",
            Self::CartnZ => "Cartn_z",
            Self::Occupancy => "occupancy",
            Self::BFactor => "B_iso_or_equiv",
            Self::FormalCharge => "pdbx_formal_charge",
            Self::ModelNum => "pdbx_PDB_model_num",
        }
    }

    /// The field an item name denotes, where the lowerer reads one.
    #[must_use]
    pub fn from_item(item: &str) -> Option<Self> {
        match item.as_bytes() {
            b"group_PDB" => Some(Self::GroupPdb),
            b"id" => Some(Self::Id),
            b"type_symbol" => Some(Self::TypeSymbol),
            b"label_atom_id" => Some(Self::LabelAtomId),
            b"auth_atom_id" => Some(Self::AuthAtomId),
            b"label_alt_id" => Some(Self::LabelAltId),
            b"label_comp_id" => Some(Self::LabelCompId),
            b"auth_comp_id" => Some(Self::AuthCompId),
            b"label_asym_id" => Some(Self::LabelAsymId),
            b"auth_asym_id" => Some(Self::AuthAsymId),
            b"label_entity_id" => Some(Self::LabelEntityId),
            b"label_seq_id" => Some(Self::LabelSeqId),
            b"auth_seq_id" => Some(Self::AuthSeqId),
            b"pdbx_PDB_ins_code" => Some(Self::InsCode),
            b"Cartn_x" => Some(Self::CartnX),
            b"Cartn_y" => Some(Self::CartnY),
            b"Cartn_z" => Some(Self::CartnZ),
            b"occupancy" => Some(Self::Occupancy),
            b"B_iso_or_equiv" => Some(Self::BFactor),
            b"pdbx_formal_charge" => Some(Self::FormalCharge),
            b"pdbx_PDB_model_num" => Some(Self::ModelNum),
            _ => None,
        }
    }

    /// Position in a row that stores fields by ordinal.
    #[must_use]
    pub const fn position(self) -> usize {
        self as usize
    }
}

/// One borrowed `atom_site` row presented to the shared structure lowerer.
///
/// Implementations must use the same value interpretation as [`crate::CifValue`]:
/// sentinels are absent, integer and floating values remain typed, and numeric
/// identifiers are rendered as their canonical decimal text.
#[doc(hidden)]
pub trait AtomSiteRow {
    /// Zero-based row number used when attaching diagnostics to the source.
    fn row(&self) -> usize;

    /// A textual value, excluding numeric values and sentinels.
    fn text(&self, field: Field) -> Option<&str>;

    /// An identifier, including numeric identifiers rendered as text.
    fn identifier(&self, field: Field) -> Option<Cow<'_, str>>;

    /// A whole-number value.
    fn integer(&self, field: Field) -> Option<i64>;

    /// A numeric value, including integers.
    fn float(&self, field: Field) -> Option<f64>;

    /// Whether the item was recorded rather than written as a sentinel.
    fn is_recorded(&self, field: Field) -> bool;
}

/// Synchronous destination for borrowed `atom_site` rows.
///
/// Rows must be fed in deposition order. Implementations consume the row and
/// all returned borrows before this method returns.
#[doc(hidden)]
pub trait AtomSiteRowSink {
    /// Lowers one complete row.
    fn feed(&mut self, row: &dyn AtomSiteRow);
}

impl AtomSiteRow for Rows<'_> {
    fn row(&self) -> usize {
        Rows::row(self)
    }

    fn text(&self, field: Field) -> Option<&str> {
        Rows::text(self, field.item())
    }

    fn identifier(&self, field: Field) -> Option<Cow<'_, str>> {
        Rows::identifier(self, field.item())
    }

    fn integer(&self, field: Field) -> Option<i64> {
        Rows::integer(self, field.item())
    }

    fn float(&self, field: Field) -> Option<f64> {
        Rows::float(self, field.item())
    }

    fn is_recorded(&self, field: Field) -> bool {
        Rows::is_recorded(self, field.item())
    }
}
