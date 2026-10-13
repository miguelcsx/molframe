//! Chain labels that name one chain each.
//!
//! A format with no separate normalised label, such as PDB, can hold several
//! chains under one identifier, a polymer and its waters after a `TER`. Formats
//! that key chains by the normalised label would merge them on the way back in,
//! so a writer asks for a view in which each chain has its own.

use super::{CoordinateStore, Structure, StructureData};
use crate::index::ChainIndex;
use crate::symbol::SymbolId;
use std::collections::HashSet;
use std::sync::Arc;

impl Structure {
    /// This structure with repeated normalised chain labels made distinct, or
    /// `None` when every chain already has its own.
    ///
    /// A repeat takes the label with the lowest unused number appended; the
    /// depositor label is untouched, so the original identifier survives in its
    /// own namespace.
    #[must_use]
    pub fn with_distinct_chain_labels(&self) -> Option<Self> {
        let mut data = self.0.as_ref().clone();
        let mut changed = distinct_in(&mut data);
        if let CoordinateStore::Ragged { models } = &mut data.coords {
            for model in models {
                if let Some(distinct) = model.with_distinct_chain_labels() {
                    *model = distinct;
                    changed = true;
                }
            }
        }
        changed.then(|| Self(Arc::new(data)))
    }
}

fn distinct_in(data: &mut StructureData) -> bool {
    let chains: Vec<(ChainIndex, SymbolId)> = data
        .topology
        .chains
        .iter()
        .filter_map(|chain| Some((chain, data.topology.chains.label_asym_id(chain)?)))
        .collect();
    let mut taken: HashSet<SymbolId> = chains.iter().map(|&(_, label)| label).collect();
    let mut seen = HashSet::new();
    let mut changed = false;
    for (chain, label) in chains {
        if seen.insert(label) {
            continue;
        }
        let base = data
            .dictionary
            .resolve(label)
            .unwrap_or_default()
            .to_owned();
        let fresh = (2u32..100_000).find_map(|number| {
            let symbol = data.dictionary.intern(&format!("{base}{number}")).ok()?;
            taken.insert(symbol).then_some(symbol)
        });
        let Some(fresh) = fresh else { continue };
        if data.topology.chains.relabel(chain, fresh).is_ok() {
            seen.insert(fresh);
            changed = true;
        }
    }
    changed
}

#[cfg(test)]
#[path = "chain_labels_tests.rs"]
mod tests;
