//! The models of a structure, in deposition order.

use crate::index::ModelIndex;
use std::ops::Range;
use std::sync::Arc;

/// The models of a structure, in deposition order.
#[derive(Clone, Debug, Default)]
pub struct ModelTable {
    first_chain: Arc<Vec<u32>>,
    chain_count: Arc<Vec<u32>>,
    model_num: Arc<Vec<i32>>,
}

impl ModelTable {
    /// The number of models.
    #[must_use]
    pub fn len(&self) -> usize {
        self.model_num.len()
    }

    /// Returns true when the structure has no models.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.model_num.is_empty()
    }

    /// Appends a model covering a range of chains.
    pub fn push(&mut self, model_num: i32, chains: Range<u32>) -> ModelIndex {
        let position = self.model_num.len() as u32;
        Arc::make_mut(&mut self.first_chain).push(chains.start);
        Arc::make_mut(&mut self.chain_count).push(chains.end.saturating_sub(chains.start));
        Arc::make_mut(&mut self.model_num).push(model_num);
        ModelIndex::new(position)
    }

    /// The number the model was deposited under.
    ///
    /// Never renumbered. A file whose models are 1, 5 and 7 keeps those numbers,
    /// because a result that cites model 5 must mean the same model the
    /// depositor called 5.
    #[must_use]
    pub fn model_num(&self, model: ModelIndex) -> Option<i32> {
        self.model_num.get(model.as_usize()).copied()
    }

    /// The chains this model contains.
    #[must_use]
    pub fn chains(&self, model: ModelIndex) -> Option<Range<u32>> {
        stored_range(&self.first_chain, &self.chain_count, model.as_usize())
    }

    /// Every model position.
    pub fn iter(&self) -> impl Iterator<Item = ModelIndex> + '_ {
        (0..self.model_num.len() as u32).map(ModelIndex::new)
    }
}

fn stored_range(starts: &[u32], counts: &[u32], index: usize) -> Option<Range<u32>> {
    let start = *starts.get(index)?;
    let count = *counts.get(index)?;

    Some(start..start.saturating_add(count))
}
