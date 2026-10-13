//! Pluggable deterministic component providers and CCD lowering.

mod bundled;
mod embedded;
mod embedded_scan;
mod items;
mod layered;
mod lazy;
mod lowering;
mod stack;
mod validation;
mod views;

use crate::model::Component;
pub use bundled::{CORE_DICTIONARY_VERSION, CoreProvider};
pub use embedded::EmbeddedProvider;
use items::one_finding;
pub use layered::LayeredProvider;
pub use lazy::LazyCifProvider;
use lowering::component;
use molframe_cif::Document;
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::InputBuffer;
pub use stack::{
    EMBEDDED_COMPONENTS_EXTENSION, embedded_components, resolution_stack, take_embedded_components,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// A source of chemical component definitions.
pub trait ComponentProvider: Send + Sync {
    /// Returns a component, or `None` when this provider does not know it.
    ///
    /// # Errors
    ///
    /// Returns a registered finding when provider storage cannot be read.
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic>;

    /// Exact dictionary version supplied by this provider.
    fn version(&self) -> &DictionaryVersion;
}

/// An in-memory deterministic provider.
#[derive(Clone, Debug)]
pub struct MemoryProvider {
    version: DictionaryVersion,
    pub(super) components: BTreeMap<Box<str>, Arc<Component>>,
}

impl MemoryProvider {
    /// Builds an index without replacing duplicate component definitions.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic when two supplied components share an identifier.
    pub fn new(
        version: DictionaryVersion,
        components: impl IntoIterator<Item = Component>,
    ) -> Result<Self, Diagnostic> {
        let mut indexed = BTreeMap::new();
        for component in components {
            let id = component.id.clone();
            if indexed.insert(id.clone(), Arc::new(component)).is_some() {
                return Err(Diagnostic::new(Code::E2005)
                    .in_category("chem_comp")
                    .with_context("id", id.to_string()));
            }
        }
        Ok(Self {
            version,
            components: indexed,
        })
    }
}

impl ComponentProvider for MemoryProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        Ok(self.components.get(component_id).cloned())
    }

    fn version(&self) -> &DictionaryVersion {
        &self.version
    }
}

/// Component provider lowered from a CCD mmCIF document.
#[derive(Clone, Debug)]
pub struct CifProvider(MemoryProvider);

impl CifProvider {
    /// Lowers every component block in a parsed CCD document.
    ///
    /// # Errors
    ///
    /// Returns every missing, malformed, or duplicate component field without
    /// constructing a partial dictionary.
    pub fn from_document(
        document: &Document,
        version: DictionaryVersion,
    ) -> Result<Self, Vec<Diagnostic>> {
        let mut components = Vec::new();
        let mut findings = Vec::new();
        let mut identifiers = std::collections::BTreeSet::new();
        for block in document.blocks() {
            match component(block) {
                Ok(component) if identifiers.insert(component.id.clone()) => {
                    components.push(component);
                }
                Ok(component) => findings.push(
                    Diagnostic::new(Code::E2005)
                        .in_category("chem_comp")
                        .with_context("id", component.id.to_string())
                        .with_context("block", block.name()),
                ),
                Err(block_findings) => findings.extend(block_findings),
            }
        }
        if findings.is_empty() {
            MemoryProvider::new(version, components)
                .map(Self)
                .map_err(one_finding)
        } else {
            Err(findings)
        }
    }
}

impl ComponentProvider for CifProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        self.0.get(component_id)
    }

    fn version(&self) -> &DictionaryVersion {
        self.0.version()
    }
}

/// Reads a CCD document and lowers every component through the shared CIF parser.
///
/// # Errors
///
/// Returns parser diagnostics or component-definition errors without producing
/// a partial provider.
pub fn read_ccd(
    input: &InputBuffer,
    version: DictionaryVersion,
) -> Result<(CifProvider, Vec<Diagnostic>), Vec<Diagnostic>> {
    let (document, findings) = molframe_cif::parse(input)?;
    CifProvider::from_document(&document, version).map(|provider| (provider, findings))
}

#[cfg(test)]
#[path = "provider_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "provider_core_tests.rs"]
mod core_tests;

#[cfg(test)]
#[path = "provider_stack_tests.rs"]
mod stack_tests;

#[cfg(test)]
#[path = "provider_embedded_tests.rs"]
mod embedded_tests;
