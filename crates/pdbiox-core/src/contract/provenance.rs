//! Where a result came from.
//!
//! A record that travels with a number and says what produced it: which version
//! of the software, which input, which policy, which reference data. It is what
//! turns "we computed 247 contacts" into something a reader can check.
//!
//! Time is deliberately excluded from the fingerprint. Two identical analyses
//! run a week apart must fingerprint identically, or the fingerprint identifies
//! the run rather than the result.

use super::policy::{AnalysisPolicy, Fingerprint, ProfileId};
use std::fmt;
use std::path::{Path, PathBuf};

/// Where an input came from.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
#[non_exhaustive]
pub enum SourceRef {
    /// Nothing was read.
    #[default]
    None,
    /// A file on disk.
    Path(PathBuf),
    /// A network location.
    Url(Box<str>),
    /// Bytes the caller supplied.
    Memory,
}

impl SourceRef {
    /// Names a file.
    #[must_use]
    pub fn path(path: impl Into<PathBuf>) -> Self {
        Self::Path(path.into())
    }

    /// The file this refers to, if it is a file.
    #[must_use]
    pub fn as_path(&self) -> Option<&Path> {
        match self {
            Self::Path(path) => Some(path),
            _ => None,
        }
    }
}

impl fmt::Display for SourceRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => f.write_str("(none)"),
            Self::Path(path) => write!(f, "{}", path.display()),
            Self::Url(url) => f.write_str(url),
            Self::Memory => f.write_str("(memory)"),
        }
    }
}

/// The version of a body of reference data.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct DictionaryVersion(Box<str>);

impl DictionaryVersion {
    /// Records a version.
    #[must_use]
    pub fn new(version: impl Into<Box<str>>) -> Self {
        Self(version.into())
    }

    /// The version string.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Everything needed to reproduce a result.
///
/// # Examples
///
/// ```
/// use pdbiox_core::contract::{AnalysisPolicy, Provenance, SourceRef};
///
/// let policy = AnalysisPolicy::default();
/// let record = Provenance::new(&policy).with_source(SourceRef::path("1abc.cif"));
///
/// assert_eq!(record.profile.map(|profile| profile.as_str()), Some("pdbiox-default-1.0"));
/// assert_eq!(record.policy_fingerprint, policy.fingerprint());
/// ```
#[derive(Clone, Debug)]
pub struct Provenance {
    /// The version of pdbiox that produced the result.
    pub pdbiox_version: &'static str,
    /// Where the input came from.
    pub input_source: SourceRef,
    /// A fingerprint of the input's bytes, where they were read.
    pub input_fingerprint: Option<Fingerprint>,
    /// The policy in force.
    pub policy: AnalysisPolicy,
    /// A fingerprint of that policy.
    pub policy_fingerprint: Fingerprint,
    /// The named profile, where the policy is one unaltered.
    pub profile: Option<ProfileId>,
    /// The version of the schema dictionary consulted.
    pub schema_version: Option<DictionaryVersion>,
    /// The version of the chemical component data consulted.
    pub component_version: Option<DictionaryVersion>,
    /// When the analysis ran, if the caller asked for it recorded.
    ///
    /// Excluded from every fingerprint.
    pub timestamp: Option<Box<str>>,
}

impl Provenance {
    /// Starts a record under `policy`.
    #[must_use]
    pub fn new(policy: &AnalysisPolicy) -> Self {
        Self {
            pdbiox_version: env!("CARGO_PKG_VERSION"),
            input_source: SourceRef::None,
            input_fingerprint: None,
            policy: policy.clone(),
            policy_fingerprint: policy.fingerprint(),
            profile: policy.profile(),
            schema_version: None,
            component_version: None,
            timestamp: None,
        }
    }

    /// Records where the input came from.
    #[must_use]
    pub fn with_source(mut self, source: SourceRef) -> Self {
        self.input_source = source;
        self
    }

    /// Records a fingerprint of the input's bytes.
    #[must_use]
    pub const fn with_input_fingerprint(mut self, fingerprint: Fingerprint) -> Self {
        self.input_fingerprint = Some(fingerprint);
        self
    }

    /// Records when the analysis ran.
    #[must_use]
    pub fn with_timestamp(mut self, timestamp: impl Into<Box<str>>) -> Self {
        self.timestamp = Some(timestamp.into());
        self
    }

    /// A fingerprint of everything that affects the result.
    ///
    /// Two runs over the same input under the same policy agree here whatever
    /// the clock said, which is what makes the value a property of the result
    /// rather than of the occasion.
    #[must_use]
    pub fn fingerprint(&self) -> Fingerprint {
        let mut bytes = Vec::with_capacity(64);
        bytes.extend_from_slice(self.pdbiox_version.as_bytes());
        bytes.extend_from_slice(&self.policy_fingerprint.get().to_le_bytes());
        if let Some(input) = self.input_fingerprint {
            bytes.extend_from_slice(&input.get().to_le_bytes());
        }
        if let Some(version) = &self.schema_version {
            bytes.extend_from_slice(version.as_str().as_bytes());
        }
        if let Some(version) = &self.component_version {
            bytes.extend_from_slice(version.as_str().as_bytes());
        }
        Fingerprint::of(&bytes)
    }
}

#[cfg(test)]
#[path = "provenance_tests.rs"]
mod tests;
