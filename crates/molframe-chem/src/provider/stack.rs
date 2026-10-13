//! The component resolution order every entry point shares.
//!
//! A structure's own definitions answer first, an explicit user dictionary
//! second and the bundled core set last, so a file that defines its ligand and
//! a user who supplies a full dictionary each win over the compact defaults
//! without any caller composing layers by hand. Definitions a file gives
//! without formal charges rank below the dictionaries, which state them.

use super::{ComponentProvider, CoreProvider, EmbeddedProvider, LayeredProvider};
use crate::model::Component;
use molframe_cif::Document;
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::structure::Structure;
use std::sync::Arc;

/// Extension key under which a reader attaches a file's own definitions.
pub const EMBEDDED_COMPONENTS_EXTENSION: &str = "molframe.chem.embedded_components.v1";

/// Version label recorded for definitions that came from the structure file.
const EMBEDDED_VERSION: &str = "embedded";

/// Lowers the definitions a parsed structure document carries, or `None` when
/// it defines no component of its own.
///
/// # Errors
///
/// Returns every malformed or duplicate embedded definition.
pub fn embedded_components(
    document: &Document,
) -> Result<Option<EmbeddedProvider>, Vec<Diagnostic>> {
    let provider =
        EmbeddedProvider::from_document(document, DictionaryVersion::new(EMBEDDED_VERSION))?;
    Ok((!provider.is_empty()).then_some(provider))
}

/// The provider for annotating `structure`: in-file definitions that state
/// charges, then `user`, then the bundled core set, then in-file definitions
/// that give connectivity only.
#[must_use]
pub fn resolution_stack(
    structure: &Structure,
    user: Option<Arc<dyn ComponentProvider>>,
) -> LayeredProvider {
    let mut layers: Vec<Arc<dyn ComponentProvider>> = Vec::new();
    let embedded = structure
        .extensions()
        .get_shared::<EmbeddedProvider>(EMBEDDED_COMPONENTS_EXTENSION);
    if let Some(embedded) = &embedded {
        layers.push(Arc::new(Tier {
            source: Arc::clone(embedded),
            charged: true,
        }));
    }
    if let Some(user) = user {
        layers.push(user);
    }
    layers.push(Arc::new(CoreProvider));
    if let Some(embedded) = embedded {
        layers.push(Arc::new(Tier {
            source: embedded,
            charged: false,
        }));
    }
    LayeredProvider::new(layers)
}

/// One half of a file's definitions: those that state charges or those that
/// do not.
struct Tier {
    source: Arc<EmbeddedProvider>,
    charged: bool,
}

impl ComponentProvider for Tier {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        if self.source.states_charges(component_id) == self.charged {
            self.source.get(component_id)
        } else {
            Ok(None)
        }
    }

    fn version(&self) -> &DictionaryVersion {
        self.source.version()
    }
}
