use super::{AnalysisPolicy, DictionaryVersion, Fingerprint, Provenance};

/// Versions of executable and reference data available for a re-execution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReexecutionEnvironment<'a> {
    /// Installed pdbiox version.
    pub pdbiox_version: &'a str,
    /// Installed schema dictionary version, when one is in use.
    pub schema_version: Option<&'a str>,
    /// Installed chemical component dictionary version, when one is in use.
    pub component_version: Option<&'a str>,
}

impl ReexecutionEnvironment<'static> {
    /// Environment for an analysis with no external dictionaries.
    #[must_use]
    pub const fn current() -> Self {
        Self {
            pdbiox_version: env!("CARGO_PKG_VERSION"),
            schema_version: None,
            component_version: None,
        }
    }
}

/// A value recomputed under an exactly validated provenance record.
#[derive(Clone, Debug)]
pub struct Reexecution<T> {
    /// Newly computed value.
    pub value: T,
    /// Exact record whose input, policy and reference versions were enforced.
    pub provenance: Provenance,
}

/// Reason exact re-execution was refused.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum ReexecutionError {
    /// The original run did not fingerprint its input bytes.
    #[error("provenance has no input fingerprint")]
    MissingInputFingerprint,
    /// Supplied bytes are not the bytes used by the original run.
    #[error("input fingerprint does not match provenance")]
    InputMismatch,
    /// The policy record was changed without updating its fingerprint.
    #[error("policy fingerprint does not match the recorded policy")]
    PolicyMismatch,
    /// Installed executable version differs from the recorded version.
    #[error("pdbiox version mismatch: expected {expected}, found {observed}")]
    SoftwareVersionMismatch {
        /// Version recorded by the original run.
        expected: Box<str>,
        /// Installed version supplied by the caller.
        observed: Box<str>,
    },
    /// Installed schema dictionary differs from the recorded version.
    #[error("schema dictionary version mismatch")]
    SchemaVersionMismatch,
    /// Installed component dictionary differs from the recorded version.
    #[error("component dictionary version mismatch")]
    ComponentVersionMismatch,
}

/// Recomputes a value only after reproducing every recorded deterministic input.
///
/// The caller supplies bytes instead of allowing hidden filesystem or network
/// access. Their fingerprint, the full policy, pdbiox version and both optional
/// dictionary versions are checked before `run` is invoked. The callback
/// receives the exact recorded policy rather than a reconstructed default.
///
/// # Errors
///
/// Refuses to invoke `run` if any recorded deterministic input is absent or
/// differs from the current environment.
pub fn reexecute_from_provenance<T>(
    provenance: &Provenance,
    input: &[u8],
    environment: ReexecutionEnvironment<'_>,
    run: impl FnOnce(&[u8], &AnalysisPolicy) -> T,
) -> Result<Reexecution<T>, ReexecutionError> {
    validate(provenance, input, environment)?;
    Ok(Reexecution {
        value: run(input, &provenance.policy),
        provenance: provenance.clone(),
    })
}

fn validate(
    provenance: &Provenance,
    input: &[u8],
    environment: ReexecutionEnvironment<'_>,
) -> Result<(), ReexecutionError> {
    let Some(expected_input) = provenance.input_fingerprint else {
        return Err(ReexecutionError::MissingInputFingerprint);
    };
    if expected_input != Fingerprint::of(input) {
        return Err(ReexecutionError::InputMismatch);
    }
    if provenance.policy_fingerprint != provenance.policy.fingerprint() {
        return Err(ReexecutionError::PolicyMismatch);
    }
    if provenance.pdbiox_version != environment.pdbiox_version {
        return Err(ReexecutionError::SoftwareVersionMismatch {
            expected: provenance.pdbiox_version.into(),
            observed: environment.pdbiox_version.into(),
        });
    }
    check_dictionary(
        provenance.schema_version.as_ref(),
        environment.schema_version,
        ReexecutionError::SchemaVersionMismatch,
    )?;
    check_dictionary(
        provenance.component_version.as_ref(),
        environment.component_version,
        ReexecutionError::ComponentVersionMismatch,
    )
}

fn check_dictionary(
    expected: Option<&DictionaryVersion>,
    observed: Option<&str>,
    mismatch: ReexecutionError,
) -> Result<(), ReexecutionError> {
    if expected.map(DictionaryVersion::as_str) == observed {
        Ok(())
    } else {
        Err(mismatch)
    }
}

#[cfg(test)]
#[path = "reexecution_tests.rs"]
mod tests;
