//! Declarative filters over manifest statistics.

use super::{DatasetError, ManifestEntry};

/// Optional predicates combined with logical AND.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DatasetFilter {
    /// Exclusive upper resolution bound in ångström.
    pub resolution_below: Option<f32>,
    /// Exact experimental method.
    pub method: Option<Box<str>>,
    /// Inclusive minimum atom count.
    pub minimum_atoms: Option<u64>,
    /// Exclusive maximum atom count.
    pub maximum_atoms: Option<u64>,
    /// Required manifest tag.
    pub tag: Option<Box<str>>,
}

impl DatasetFilter {
    pub(crate) fn validate(&self) -> Result<(), DatasetError> {
        if self
            .resolution_below
            .is_some_and(|value| !value.is_finite() || value <= 0.0)
            || self
                .minimum_atoms
                .zip(self.maximum_atoms)
                .is_some_and(|(minimum, maximum)| minimum >= maximum)
        {
            Err(DatasetError::InvalidFilter)
        } else {
            Ok(())
        }
    }

    pub(crate) fn matches(&self, entry: &ManifestEntry) -> bool {
        self.resolution_below
            .is_none_or(|limit| entry.resolution.is_some_and(|value| value < limit))
            && self
                .method
                .as_deref()
                .is_none_or(|method| entry.method.as_deref() == Some(method))
            && self
                .minimum_atoms
                .is_none_or(|minimum| entry.atom_count >= minimum)
            && self
                .maximum_atoms
                .is_none_or(|maximum| entry.atom_count < maximum)
            && self
                .tag
                .as_deref()
                .is_none_or(|tag| entry.tags.iter().any(|value| value.as_ref() == tag))
    }
}
