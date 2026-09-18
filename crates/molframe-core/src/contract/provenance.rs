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
use std::collections::BTreeMap;
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

/// Stable name and version of the algorithm that produced a result.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct AlgorithmId {
    name: Box<str>,
    version: Box<str>,
}

impl AlgorithmId {
    /// Creates a versioned algorithm identifier.
    #[must_use]
    pub fn new(name: impl Into<Box<str>>, version: impl Into<Box<str>>) -> Self {
        Self {
            name: name.into(),
            version: version.into(),
        }
    }

    /// Stable algorithm name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Algorithm version, independent of the crate release.
    #[must_use]
    pub fn version(&self) -> &str {
        &self.version
    }
}

/// A typed, deterministic result-affecting algorithm parameter.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum ParameterValue {
    /// Boolean control.
    Boolean(bool),
    /// Signed integer control.
    Integer(i64),
    /// IEEE-754 value stored by its exact bits.
    Float(u64),
    /// Textual choice or identifier.
    Text(Box<str>),
}

impl ParameterValue {
    /// Records an exact finite floating-point parameter.
    #[must_use]
    pub fn finite_float(value: f64) -> Option<Self> {
        value.is_finite().then_some(Self::Float(value.to_bits()))
    }

    pub(crate) fn serialised(&self) -> String {
        match self {
            Self::Boolean(value) => format!("bool:{value}"),
            Self::Integer(value) => format!("int:{value}"),
            Self::Float(bits) => format!("float:{bits:016x}"),
            Self::Text(value) => format!("text:{value}"),
        }
    }
}

/// Sorted algorithm parameters that participate in provenance fingerprints.
pub type AnalysisParameters = BTreeMap<Box<str>, ParameterValue>;

/// Everything needed to reproduce a result.
///
/// # Examples
///
/// ```
/// use molframe_core::contract::{AnalysisPolicy, Provenance, SourceRef};
///
/// let policy = AnalysisPolicy::default();
/// let record = Provenance::new(&policy).with_source(SourceRef::path("1abc.cif"));
///
/// assert_eq!(record.profile.map(|profile| profile.as_str()), Some("molframe-default-1.0"));
/// assert_eq!(record.policy_fingerprint, policy.fingerprint());
/// ```
#[derive(Clone, Debug)]
pub struct Provenance {
    /// The version of molframe that produced the result.
    pub molframe_version: &'static str,
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
    /// Versioned algorithm, when the result is algorithm-specific.
    pub algorithm: Option<AlgorithmId>,
    /// Sorted typed parameters that can change the result.
    pub parameters: AnalysisParameters,
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
            molframe_version: env!("CARGO_PKG_VERSION"),
            input_source: SourceRef::None,
            input_fingerprint: None,
            policy: policy.clone(),
            policy_fingerprint: policy.fingerprint(),
            profile: policy.profile(),
            algorithm: None,
            parameters: BTreeMap::new(),
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

    /// Records the versioned algorithm used by this analysis.
    #[must_use]
    pub fn with_algorithm(mut self, algorithm: AlgorithmId) -> Self {
        self.algorithm = Some(algorithm);
        self
    }

    /// Records or replaces one result-affecting parameter.
    #[must_use]
    pub fn with_parameter(mut self, name: impl Into<Box<str>>, value: ParameterValue) -> Self {
        self.parameters.insert(name.into(), value);
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
        bytes.extend_from_slice(self.molframe_version.as_bytes());
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
        if let Some(algorithm) = &self.algorithm {
            fingerprint_field(&mut bytes, b"algorithm.name", algorithm.name());
            fingerprint_field(&mut bytes, b"algorithm.version", algorithm.version());
        }
        for (name, value) in &self.parameters {
            fingerprint_field(&mut bytes, name.as_bytes(), &value.serialised());
        }
        Fingerprint::of(&bytes)
    }
}

fn fingerprint_field(bytes: &mut Vec<u8>, name: &[u8], value: &str) {
    bytes.extend_from_slice(&(name.len() as u64).to_le_bytes());
    bytes.extend_from_slice(name);
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

#[cfg(test)]
#[path = "provenance_tests.rs"]
mod tests;
