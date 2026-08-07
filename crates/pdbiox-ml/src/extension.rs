//! Domain extension metadata over standard Arrow storage types.

use arrow::datatypes::{DataType, Field};
use std::collections::HashMap;
use std::sync::Arc;

const EXTENSION_NAME_KEY: &str = "ARROW:extension:name";
const EXPORT_COST_KEY: &str = "pdbiox:export_cost";

/// Whether an exported Arrow column aliases storage or materialises it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExportCost {
    /// Arrow sees the immutable pdbiox allocation directly.
    ZeroCopy,
    /// An encoded or derived column was decoded in one pass.
    Decode,
}

/// Domain extension types published by the Arrow adapter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PdbioxExtension {
    /// Index into the atom table.
    AtomIndex,
    /// Index into the residue table.
    ResidueIndex,
    /// Index into the chain table.
    ChainIndex,
    /// Index into the entity table.
    EntityIndex,
    /// Cartesian coordinate triple in angstroms.
    Coordinates3f,
    /// Atomic number.
    Element,
    /// Identifier in the structure symbol dictionary.
    SymbolId,
    /// Alternate-location symbol identifier; zero is blank.
    Altloc,
    /// Serialised adaptive atom selection.
    Selection,
    /// Boolean companion distinguishing unknown from inapplicable nulls.
    Validity,
}

impl PdbioxExtension {
    /// Stable extension name stored in Arrow field metadata.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::AtomIndex => "pdbiox.atom_index",
            Self::ResidueIndex => "pdbiox.residue_index",
            Self::ChainIndex => "pdbiox.chain_index",
            Self::EntityIndex => "pdbiox.entity_index",
            Self::Coordinates3f => "pdbiox.coordinates3f",
            Self::Element => "pdbiox.element",
            Self::SymbolId => "pdbiox.symbol_id",
            Self::Altloc => "pdbiox.altloc",
            Self::Selection => "pdbiox.selection",
            Self::Validity => "pdbiox.validity",
        }
    }

    /// Standard Arrow storage type visible to unaware consumers.
    #[must_use]
    pub fn storage_type(self) -> DataType {
        match self {
            Self::AtomIndex
            | Self::ResidueIndex
            | Self::ChainIndex
            | Self::EntityIndex
            | Self::SymbolId
            | Self::Altloc => DataType::UInt32,
            Self::Coordinates3f => {
                DataType::FixedSizeList(Arc::new(Field::new("item", DataType::Float32, false)), 3)
            }
            Self::Element => DataType::UInt8,
            Self::Selection => DataType::Binary,
            Self::Validity => DataType::Boolean,
        }
    }
}

impl ExportCost {
    const fn label(self) -> &'static str {
        match self {
            Self::ZeroCopy => "zero-copy",
            Self::Decode => "decode",
        }
    }
}

pub(crate) fn field(
    name: &str,
    data_type: DataType,
    nullable: bool,
    extension: Option<&str>,
    cost: ExportCost,
) -> Field {
    let mut metadata = HashMap::from([(EXPORT_COST_KEY.to_owned(), cost.label().to_owned())]);
    if let Some(extension) = extension {
        metadata.insert(EXTENSION_NAME_KEY.to_owned(), extension.to_owned());
    }
    Field::new(name, data_type, nullable).with_metadata(metadata)
}

/// Returns the pdbiox extension name carried by a field.
#[must_use]
pub fn extension_name(field: &Field) -> Option<&str> {
    field.metadata().get(EXTENSION_NAME_KEY).map(String::as_str)
}

#[cfg(test)]
#[path = "extension_tests.rs"]
mod tests;
