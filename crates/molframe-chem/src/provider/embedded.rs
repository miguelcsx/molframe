//! Component definitions carried inside a structure file.
//!
//! A structure mmCIF may embed `chem_comp_atom` and `chem_comp_bond` rows for
//! the components it uses, interleaved in one data block and keyed by
//! `comp_id`. This lowers those rows into a provider. A component listed only
//! in `chem_comp`, with no atom rows, is not defined: a ligand with no atoms
//! would silently read as complete chemistry, so it stays unknown.

use super::lowering::{Rows, component_from};
use super::{ComponentProvider, MemoryProvider, one_finding};
use crate::Component;
use molframe_cif::{Category, DataBlock, Document};
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::Diagnostic;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Definitions of the components a structure file defines for itself.
#[derive(Clone, Debug)]
pub struct EmbeddedProvider(MemoryProvider);

impl EmbeddedProvider {
    /// Lowers every component with embedded atom rows in a parsed structure
    /// document.
    ///
    /// The result may be empty, which is the usual case for a deposited file
    /// that lists components without defining them; see [`Self::is_empty`].
    ///
    /// # Errors
    ///
    /// Returns every malformed or duplicate embedded definition without
    /// constructing a partial provider.
    pub fn from_document(
        document: &Document,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let mut components = Vec::new();
        let mut findings = Vec::new();
        for block in document.blocks() {
            lower_block(block, &mut components, &mut findings);
        }
        if !findings.is_empty() {
            return Err(findings);
        }
        MemoryProvider::new(version, components)
            .map(Self)
            .map_err(one_finding)
    }

    /// Whether the file defines no component of its own.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.components.is_empty()
    }
}

fn lower_block(block: &DataBlock, components: &mut Vec<Component>, findings: &mut Vec<Diagnostic>) {
    let (Some(metadata), Some(atoms)) = (
        block.category("chem_comp"),
        block.category("chem_comp_atom"),
    ) else {
        return;
    };
    let atom_rows = group_rows(atoms);
    let bond_rows = block.category("chem_comp_bond").map(group_rows);
    for row in 0..metadata.row_count() {
        let Some(id) = metadata.identifier("id", row) else {
            continue;
        };
        let Some(rows) = atom_rows.get(id.as_ref()) else {
            continue;
        };
        let bonds = bond_rows.as_ref().and_then(|groups| {
            let rows = groups.get(id.as_ref())?;
            block
                .category("chem_comp_bond")
                .map(|category| Rows { category, rows })
        });
        match component_from(
            metadata,
            row,
            Some(Rows {
                category: atoms,
                rows,
            }),
            bonds,
        ) {
            Ok(component) => components.push(component),
            Err(found) => findings.extend(found),
        }
    }
}

fn group_rows(category: &Category) -> BTreeMap<String, Vec<usize>> {
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for row in 0..category.row_count() {
        if let Some(id) = category.identifier("comp_id", row) {
            groups.entry(id.into_owned()).or_default().push(row);
        }
    }
    groups
}

impl ComponentProvider for EmbeddedProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        self.0.get(component_id)
    }

    fn version(&self) -> &DictionaryVersion {
        self.0.version()
    }
}
