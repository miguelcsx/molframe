//! Dense ensembles as one `(models, atoms, 3)` block, and the refusal ragged ones get.
//!
//! A dense ensemble keeps one coordinate block per model over a single topology,
//! so asking for all of them is a gather into one exactly sized buffer. A ragged
//! ensemble has no shared atom axis: there is nothing to gather into, and the
//! honest answer is a typed refusal that points at the per-model structures.

use super::Structure;
use crate::diagnostic::{Code, Diagnostic};
use crate::index::ModelIndex;

/// The positions of several models of one dense ensemble, model-major.
#[derive(Clone, Debug, PartialEq)]
pub struct EnsembleBlock {
    models: usize,
    atoms: usize,
    coordinates: Vec<f32>,
}

impl EnsembleBlock {
    /// Models in the block.
    #[must_use]
    pub const fn models(&self) -> usize {
        self.models
    }

    /// Atoms in every model.
    #[must_use]
    pub const fn atoms(&self) -> usize {
        self.atoms
    }

    /// `models * atoms * 3` values: model, then atom, then axis.
    #[must_use]
    pub fn as_flat(&self) -> &[f32] {
        &self.coordinates
    }

    /// Hands the buffer over without copying it.
    #[must_use]
    pub fn into_flat(self) -> Vec<f32> {
        self.coordinates
    }

    /// One model's positions.
    #[must_use]
    pub fn model(&self, index: usize) -> Option<&[[f32; 3]]> {
        let width = self.atoms.checked_mul(3)?;
        let start = index.checked_mul(width)?;
        let flat = self.coordinates.get(start..start.checked_add(width)?)?;
        Some(flat.as_chunks::<3>().0)
    }
}

/// The message every refusal of a ragged ensemble carries.
const RAGGED_MESSAGE: &str = "ragged ensemble: its models hold different atoms and have no shared atom axis; \
     iterate `ragged_models` (Python `structure.ragged_models`, Rust `Structure::ragged_models()`) \
     and analyse each model on its own";

impl Structure {
    /// The refusal a ragged ensemble gets from every operation that needs one atom axis.
    pub(crate) fn ragged_refusal(&self) -> Diagnostic {
        Diagnostic::new(Code::E6003)
            .with_message(RAGGED_MESSAGE)
            .with_context("reason", RAGGED_MESSAGE)
            .with_context("models", self.model_count().to_string())
    }

    /// Gathers the positions of `models` (every model when `None`) into one buffer.
    ///
    /// The buffer is allocated once at its final size, so peak memory is the
    /// result itself. Models keep the order given.
    ///
    /// # Errors
    ///
    /// Returns [`Code::E6003`] for a ragged ensemble (naming `ragged_models`) or
    /// for a model index the structure does not have, and [`Code::E7001`] when
    /// the buffer size overflows `usize` or cannot be reserved.
    pub fn ensemble_positions(
        &self,
        models: Option<&[ModelIndex]>,
    ) -> Result<EnsembleBlock, Diagnostic> {
        self.require_resolved_model()?;
        let count = self.model_count();
        let atoms = self.positions().len();
        let all: Vec<ModelIndex> = if models.is_none() {
            (0..count)
                .map(|index| {
                    ModelIndex::new(match u32::try_from(index) {
                        Ok(index) => index,
                        Err(_) => u32::MAX,
                    })
                })
                .collect()
        } else {
            Vec::new()
        };
        let chosen = match models {
            Some(models) => models,
            None => &all,
        };
        let too_big = || {
            Diagnostic::new(Code::E7001)
                .with_message("ensemble coordinates exceed addressable memory")
                .with_context("models", chosen.len().to_string())
                .with_context("atoms", atoms.to_string())
        };
        let total = chosen
            .len()
            .checked_mul(atoms)
            .and_then(|rows| rows.checked_mul(3))
            .ok_or_else(too_big)?;
        let mut coordinates = Vec::new();
        coordinates
            .try_reserve_exact(total)
            .map_err(|_| too_big())?;
        for model in chosen {
            let Some(positions) = self.model_positions(*model) else {
                return Err(Diagnostic::new(Code::E6003)
                    .with_message(format!(
                        "model {} does not exist: the structure has {count} models, addressed by position from 0",
                        model.get()
                    ))
                    .with_context("model", model.get().to_string())
                    .with_context("models", count.to_string()));
            };
            coordinates.extend(positions.iter().flatten());
        }
        Ok(EnsembleBlock {
            models: chosen.len(),
            atoms,
            coordinates,
        })
    }
}
