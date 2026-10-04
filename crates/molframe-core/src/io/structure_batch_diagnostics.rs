//! Stable diagnostic codes for the bounded batch readers' error.

use super::StructureBatchError;
use crate::diagnostic::{Code, Diagnostic};

impl From<&StructureBatchError> for Diagnostic {
    /// A wrapped diagnostic is returned as it is, so its code, remedy and span
    /// survive; the rest name the kind of limit that stopped the read.
    fn from(error: &StructureBatchError) -> Self {
        let code = match error {
            StructureBatchError::Diagnostic(inner) => return inner.clone(),
            StructureBatchError::DemandTooSmall { .. } => Code::E5101,
            StructureBatchError::RecordExceedsBudget { .. } | StructureBatchError::Memory(_) => {
                Code::E7001
            }
            StructureBatchError::Identity(_) | StructureBatchError::DictionaryFull => Code::E1903,
        };
        Self::new(code).with_message(error.to_string())
    }
}

impl From<StructureBatchError> for Diagnostic {
    fn from(error: StructureBatchError) -> Self {
        Self::from(&error)
    }
}

#[cfg(test)]
#[path = "structure_batch_diagnostics_tests.rs"]
mod tests;
