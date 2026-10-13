//! A component dictionary chosen through the table store.
//!
//! The cache holds dictionaries the user fetched (`<name>.cif`); the profile's
//! `components:` line says which one answers when the caller names none. The
//! entry's file stem is the dictionary version recorded in results. Nothing is
//! compiled in, and an entry the profile names but the cache lacks is a typed
//! [`TableError`], never a quiet fall-through.
//!
//! Cost: indexing is one linear scan of the file when first selected; the store
//! memoises it, so later selections are a `stat` and a hash probe.

use super::{ComponentProvider, LazyCifProvider};
use crate::model::Component;
use crate::{Provided, TableError, TableIdentity, TableKind, TableStore, TableTextError};
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::io::InputBuffer;
use std::path::Path;
use std::sync::Arc;

/// Table kind of a Chemical Component Dictionary file.
#[derive(Debug)]
pub struct ComponentDictionary {
    provider: LazyCifProvider,
}

impl TableKind for ComponentDictionary {
    const KIND: &'static str = "components";
    const EXTENSION: &'static str = "cif";

    fn decode(bytes: &[u8], fallback: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let stem = match Path::new(fallback.source.as_ref()).file_stem() {
            Some(stem) => stem.to_string_lossy().into_owned(),
            None => fallback.version.to_string(),
        };
        let identity = TableIdentity::supplied(&stem, &stem, &fallback.source, &fallback.licence);
        let provider = LazyCifProvider::new(
            InputBuffer::from_bytes(bytes.to_vec()),
            DictionaryVersion::new(stem),
        );
        Ok(Provided {
            identity,
            table: Self { provider },
        })
    }
}

/// A selected dictionary as a component layer.
struct Selected(Arc<Provided<ComponentDictionary>>);

impl ComponentProvider for Selected {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        self.0.table.provider.get(component_id)
    }

    fn version(&self) -> &DictionaryVersion {
        self.0.table.provider.version()
    }
}

/// The dictionary the profile selects, or `None` when the profile names none.
///
/// # Errors
///
/// [`TableError`] when the profile names an entry the cache lacks or that
/// cannot be read.
pub fn profile_dictionary(
    store: &TableStore,
) -> Result<Option<Arc<dyn ComponentProvider>>, TableError> {
    Ok(store
        .select_optional::<ComponentDictionary>("--ccd")?
        .map(|chosen| Arc::new(Selected(chosen.provided)) as Arc<dyn ComponentProvider>))
}
