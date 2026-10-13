//! Thread-safe, once-per-id access to component dictionary entries.
//!
//! A dictionary lookup parses and indexes a component, so each id is resolved
//! once and shared across every residue and worker thread.

use super::component_index::ComponentIndex;
use crate::{Component, ComponentProvider};
use molframe_core::diagnostic::Diagnostic;
use std::collections::{BTreeMap, HashMap};
use std::sync::{Arc, Mutex};

/// A component with its name index and hydrogen-bonding roles.
///
/// Both are pure functions of the component, so they are built once per
/// component id and shared by every residue and every thread.
pub(super) struct Resolved {
    pub(super) component: Arc<Component>,
    pub(super) index: ComponentIndex,
}

/// Thread-safe, once-per-id front for the component provider.
pub(super) struct Resolver<'p> {
    provider: &'p dyn ComponentProvider,
    shared: Mutex<HashMap<Box<str>, Option<Arc<Resolved>>>>,
}

impl<'p> Resolver<'p> {
    pub(super) fn new(provider: &'p dyn ComponentProvider) -> Self {
        Self {
            provider,
            shared: Mutex::new(HashMap::new()),
        }
    }

    pub(super) fn lookup(&self, id: &str) -> Result<Option<Arc<Resolved>>, Diagnostic> {
        let cached = match self.shared.lock() {
            Ok(guard) => guard.get(id).cloned(),
            Err(poisoned) => poisoned.into_inner().get(id).cloned(),
        };
        if let Some(hit) = cached {
            return Ok(hit);
        }
        // The provider runs outside the lock so distinct components lower in
        // parallel. A rare duplicate lookup is harmless: the first insert wins.
        let resolved = self.provider.get(id)?.map(|component| {
            Arc::new(Resolved {
                index: ComponentIndex::build(Arc::clone(&component)),
                component,
            })
        });
        let mut shared = match self.shared.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        };
        Ok(shared.entry(id.into()).or_insert(resolved).clone())
    }

    /// Every component resolved so far, by id.
    pub(super) fn components(&self) -> BTreeMap<Box<str>, Arc<Component>> {
        let shared = match self.shared.lock() {
            Ok(guard) => guard,
            Err(poisoned) => poisoned.into_inner(),
        };
        shared
            .iter()
            .filter_map(|(id, resolved)| {
                resolved
                    .as_ref()
                    .map(|resolved| (id.clone(), Arc::clone(&resolved.component)))
            })
            .collect()
    }
}
