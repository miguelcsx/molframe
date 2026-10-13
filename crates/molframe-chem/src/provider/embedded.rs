//! Component definitions carried inside a structure file.
//!
//! A structure mmCIF may embed `chem_comp_atom` and `chem_comp_bond` rows for
//! the components it uses, interleaved in one data block and keyed by
//! `comp_id`. The provider keeps those tables as read and lowers a component
//! only when something asks for it, then caches it. A component listed only in
//! `chem_comp`, with no atom rows, is not defined: a ligand with no atoms
//! would silently read as complete chemistry, so it stays unknown.
//!
//! Construction is `O(rows)` with no per-row allocation: one pass indexes the
//! `comp_id` columns and one borrowed pass checks every row. A malformed or
//! duplicate definition therefore still fails construction with the same
//! findings as lowering everything, so a reader reports it at read time and a
//! later lookup cannot fail. A lookup costs `O(log k)` plus, on the first ask
//! for a component, the lowering of that component's own rows.

use super::embedded_scan::Source;
use super::{ComponentProvider, one_finding};
use crate::Component;
use molframe_cif::{Category, Document};
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Code, Diagnostic};
use std::collections::BTreeSet;
use std::sync::Arc;

/// Definitions of the components a structure file defines for itself.
#[derive(Clone, Debug)]
pub struct EmbeddedProvider {
    version: DictionaryVersion,
    sources: Arc<[Source]>,
}

impl EmbeddedProvider {
    /// Indexes every component with embedded atom rows in a parsed structure
    /// document, copying the component tables out of it.
    ///
    /// The result may be empty, which is the usual case for a deposited file
    /// that lists components without defining them; see [`Self::is_empty`].
    /// A reader that owns the document should prefer [`Self::take_from`].
    ///
    /// # Errors
    ///
    /// Returns every malformed or duplicate embedded definition without
    /// constructing a partial provider.
    pub fn from_document(
        document: &Document,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let tables = document.blocks().filter_map(|block| {
            Some((
                block.category("chem_comp")?.clone(),
                block.category("chem_comp_atom")?.clone(),
                block.category("chem_comp_bond").cloned(),
            ))
        });
        Self::from_tables(tables, version)
    }

    /// Like [`Self::from_document`], moving the component tables out of the
    /// document instead of copying them. The document no longer carries them
    /// afterwards.
    ///
    /// # Errors
    ///
    /// Returns every malformed or duplicate embedded definition without
    /// constructing a partial provider.
    pub fn take_from(
        document: &mut Document,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let tables: Vec<_> = document
            .blocks_mut()
            .filter_map(|block| {
                let metadata = block.take_category("chem_comp")?;
                let atoms = block.take_category("chem_comp_atom")?;
                Some((metadata, atoms, block.take_category("chem_comp_bond")))
            })
            .collect();
        Self::from_tables(tables, version)
    }

    fn from_tables(
        tables: impl IntoIterator<Item = (Category, Category, Option<Category>)>,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let sources: Arc<[Source]> = tables
            .into_iter()
            .filter_map(|(metadata, atoms, bonds)| Source::new(metadata, atoms, bonds))
            .collect();
        if sources.iter().all(Source::is_valid) && !has_duplicate(&sources) {
            return Ok(Self { version, sources });
        }
        let mut findings: Vec<Diagnostic> = sources.iter().flat_map(Source::findings).collect();
        if findings.is_empty() {
            let mut seen = BTreeSet::new();
            let duplicate = sources
                .iter()
                .find_map(|source| source.first_duplicate(&mut seen));
            findings = one_finding(
                Diagnostic::new(Code::E2005)
                    .in_category("chem_comp")
                    .with_context("id", duplicate.map(|id| id.to_string()).unwrap_or_default()),
            );
        }
        Err(findings)
    }

    /// Whether the file stated this component's formal charges.
    ///
    /// Deposited entries list every component's atoms and bonds but leave the
    /// charges out, so those definitions describe connectivity only and a
    /// dictionary that does state charges is the better source for them.
    #[must_use]
    pub fn states_charges(&self, component_id: &str) -> bool {
        self.sources
            .iter()
            .find(|source| source.entry(component_id).is_some())
            .is_none_or(|source| source.charged)
    }

    /// Whether the file defines no component of its own.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sources.is_empty()
    }
}

/// Whether two defined components share an id, within or across blocks.
fn has_duplicate(sources: &[Source]) -> bool {
    let within = sources.iter().any(|source| {
        source
            .entries
            .windows(2)
            .any(|pair| pair[0].id == pair[1].id)
    });
    within
        || (sources.len() > 1 && {
            let mut seen = BTreeSet::new();
            sources
                .iter()
                .flat_map(|source| source.entries.iter())
                .any(|entry| !seen.insert(entry.id.as_ref()))
        })
}

impl ComponentProvider for EmbeddedProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        for source in self.sources.iter() {
            if let Some(entry) = source.entry(component_id) {
                return entry.component(source).map(Some);
            }
        }
        Ok(None)
    }

    fn version(&self) -> &DictionaryVersion {
        &self.version
    }
}
