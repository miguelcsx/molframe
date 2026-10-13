//! The compact component set shipped inside the crate.
//!
//! Covers the standard and common modified residues, nucleotides, water,
//! common ions, crystallisation additives, frequent cofactors and the
//! carbohydrates molframe names. Everything else is unknown here and must
//! come from a full dictionary or from the structure file. Definitions carry
//! atoms, charges, aromaticity and bond orders; ideal and model coordinates
//! are not bundled, so they read as absent.

use super::{ComponentProvider, LazyCifProvider};
use crate::model::Component;
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::io::InputBuffer;
use std::sync::{Arc, OnceLock};

const CORE_CIF: &[u8] = include_bytes!("../../data/core-components.cif");

/// Version label of the bundled set: the wwPDB release it was cut from.
pub const CORE_DICTIONARY_VERSION: &str = "wwpdb-ccd-2026-10-03-core";

/// The bundled core component dictionary.
///
/// Blocks are indexed once per process and lowered on first request.
#[derive(Clone, Copy, Debug, Default)]
pub struct CoreProvider;

fn shared() -> &'static LazyCifProvider {
    static SHARED: OnceLock<LazyCifProvider> = OnceLock::new();
    SHARED.get_or_init(|| {
        LazyCifProvider::new(
            InputBuffer::from_bytes(CORE_CIF.to_vec()),
            DictionaryVersion::new(CORE_DICTIONARY_VERSION),
        )
    })
}

impl CoreProvider {
    /// Number of components in the bundled set.
    #[must_use]
    pub fn component_count(&self) -> usize {
        shared().block_count()
    }
}

impl ComponentProvider for CoreProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        shared().get(component_id)
    }

    fn version(&self) -> &DictionaryVersion {
        shared().version()
    }
}
