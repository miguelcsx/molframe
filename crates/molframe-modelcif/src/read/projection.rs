//! Reusable compact projection over synchronous CIF scalar callbacks.

use super::category::ModelBuilder;
use crate::{ModelCif, ModelCifError, ModelCifOptions};
use molframe_cif::{CifEventSink, CifScalar};
use molframe_core::{Diagnostic, span::ByteSpan};

/// Event projection used to combine structure and `ModelCIF` reading in one lex.
#[doc(hidden)]
pub struct ModelCifProjection {
    builder: ModelBuilder,
    blocks: usize,
    limit: usize,
    error: Option<ModelCifError>,
}

impl std::fmt::Debug for ModelCifProjection {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ModelCifProjection")
            .field("blocks", &self.blocks)
            .field("limit", &self.limit)
            .field("error", &self.error)
            .finish_non_exhaustive()
    }
}

impl ModelCifProjection {
    /// Creates a direct projection after validating its resource policy.
    ///
    /// # Errors
    ///
    /// Returns an invalid-limit error before consuming input.
    pub fn new(options: ModelCifOptions) -> Result<Self, ModelCifError> {
        let options = options.validate()?;
        Ok(Self {
            builder: ModelBuilder::new(options.memory_limit_bytes),
            blocks: 0,
            limit: options.memory_limit_bytes,
            error: None,
        })
    }

    fn record(&mut self, result: Result<(), ModelCifError>) {
        if self.error.is_none()
            && let Err(error) = result
        {
            self.builder = ModelBuilder::new(self.limit);
            self.error = Some(error);
        }
    }
}

impl CifEventSink for ModelCifProjection {
    type Output = Result<(ModelCif, Vec<Diagnostic>), ModelCifError>;

    fn block(&mut self, _name: &str) {
        self.blocks = self.blocks.saturating_add(1);
    }

    fn accepts_category(&self, category: &str) -> bool {
        self.blocks == 1 && category.starts_with("ma_") && self.error.is_none()
    }

    fn value(&mut self, category: &str, item: &str, value: CifScalar<'_>, _span: ByteSpan) {
        if self.blocks != 1 || !category.starts_with("ma_") || self.error.is_some() {
            return;
        }
        let result = self.builder.push(category, item, value);
        self.record(result);
    }

    fn finish(self) -> Self::Output {
        if let Some(error) = self.error {
            return Err(error);
        }
        let model = self.builder.finish()?;
        let findings = crate::lower::validate_model(&model);
        Ok((model, findings))
    }
}
