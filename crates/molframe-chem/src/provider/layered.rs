//! An ordered stack of component providers.

use super::ComponentProvider;
use crate::model::Component;
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::Diagnostic;
use std::sync::Arc;

/// Providers consulted in order; the first that knows a component answers.
///
/// A component is never merged across layers: its atoms and bonds come whole
/// from one dictionary. The version string names every layer in order, so a
/// result records that more than one source was consulted.
pub struct LayeredProvider {
    layers: Vec<Arc<dyn ComponentProvider>>,
    version: DictionaryVersion,
}

impl std::fmt::Debug for LayeredProvider {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("LayeredProvider")
            .field("version", &self.version)
            .field("layers", &self.layers.len())
            .finish()
    }
}

impl LayeredProvider {
    /// Stacks `layers`, highest priority first.
    #[must_use]
    pub fn new(layers: Vec<Arc<dyn ComponentProvider>>) -> Self {
        let names: Vec<&str> = layers
            .iter()
            .map(|layer| layer.version().as_str())
            .collect();
        let version = DictionaryVersion::new(names.join("+"));
        Self { layers, version }
    }
}

impl ComponentProvider for LayeredProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        for layer in &self.layers {
            if let Some(component) = layer.get(component_id)? {
                return Ok(Some(component));
            }
        }
        Ok(None)
    }

    fn version(&self) -> &DictionaryVersion {
        &self.version
    }
}
