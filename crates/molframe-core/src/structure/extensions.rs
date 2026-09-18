//! Type-safe, structure-local domain extensions.

use std::any::Any;
use std::collections::BTreeMap;
use std::fmt;
use std::panic::RefUnwindSafe;
use std::sync::Arc;

trait ExtensionValue: Send + Sync + RefUnwindSafe {
    fn as_any(&self) -> &dyn Any;
    fn into_any(self: Arc<Self>) -> Arc<dyn Any + Send + Sync>;
}

impl<T> ExtensionValue for T
where
    T: Any + Send + Sync + RefUnwindSafe,
{
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn into_any(self: Arc<Self>) -> Arc<dyn Any + Send + Sync> {
        self
    }
}

type SharedExtension = Arc<dyn ExtensionValue>;

/// Cold, typed metadata owned by a structure snapshot.
///
/// Keys are stable schema identifiers. Values are immutable and reference
/// counted, so cloning a structure does not copy domain metadata.
#[derive(Clone, Default)]
pub struct ExtensionStore {
    values: BTreeMap<&'static str, SharedExtension>,
}

impl ExtensionStore {
    /// Creates an empty store.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: BTreeMap::new(),
        }
    }

    /// Inserts or replaces a typed extension.
    pub fn insert<T>(&mut self, key: &'static str, value: T)
    where
        T: Any + Send + Sync + RefUnwindSafe,
    {
        let _ = self.values.insert(key, Arc::new(value));
    }

    /// Borrows an extension when its key and concrete type both match.
    #[must_use]
    pub fn get<T>(&self, key: &str) -> Option<&T>
    where
        T: Any + Send + Sync + RefUnwindSafe,
    {
        self.values.get(key)?.as_ref().as_any().downcast_ref()
    }

    /// Clones the shared handle to an extension when its type matches.
    #[must_use]
    pub fn get_shared<T>(&self, key: &str) -> Option<Arc<T>>
    where
        T: Any + Send + Sync + RefUnwindSafe,
    {
        Arc::clone(self.values.get(key)?).into_any().downcast().ok()
    }

    /// Removes every extension.
    pub fn clear(&mut self) {
        self.values.clear();
    }

    /// Returns whether no extensions are attached.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Returns the number of attached extensions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Iterates stable schema keys in lexical order.
    pub fn keys(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.values.keys().copied()
    }
}

impl fmt::Debug for ExtensionStore {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ExtensionStore")
            .field("keys", &self.values.keys().collect::<Vec<_>>())
            .finish()
    }
}
