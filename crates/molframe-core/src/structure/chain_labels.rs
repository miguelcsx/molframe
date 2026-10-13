//! Label identities that name one chain and one residue each.
//!
//! A format with no separate normalised label, such as PDB, can hold several
//! chains under one identifier, a polymer and its waters after a `TER`. Formats
//! that key chains by the normalised label would merge them on the way back in,
//! so a writer asks for a view in which each chain has its own. Insertion
//! codes raise the same problem for residues: `10` and `10A` share a number in
//! a format whose normalised sequence position carries no code.

use super::{CoordinateStore, Structure, StructureData};
use crate::index::ChainIndex;
use crate::symbol::SymbolId;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

impl Structure {
    /// This structure with repeated normalised chain labels and repeated
    /// normalised residue numbers made distinct, or `None` when every chain and
    /// residue already has its own.
    ///
    /// A repeated chain takes the label with the lowest unused number
    /// appended. Within a chain that repeats a residue number, each later
    /// residue steps the numbering on by one so the order is kept. The
    /// depositor identifiers are untouched, so the original survives in its own
    /// namespace.
    #[must_use]
    pub fn with_distinct_labels(&self) -> Option<Self> {
        let mut data = self.0.as_ref().clone();
        let mut changed = distinct_in(&mut data);
        changed |= distinct_residues(&mut data);
        if let CoordinateStore::Ragged { models } = &mut data.coords {
            for model in models {
                if let Some(distinct) = model.with_distinct_labels() {
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
            .map_or_else(String::new, str::to_owned);
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

fn distinct_residues(data: &mut StructureData) -> bool {
    let mut renumbered = Vec::new();
    for chain in data.chains() {
        let numbers: Vec<_> = chain
            .residues()
            .filter_map(|residue| {
                let code = residue.ins_code().map_or_else(String::new, str::to_owned);
                Some((residue.index(), residue.label_seq_id()?, code))
            })
            .collect();
        // Microheterogeneity repeats a number under one insertion code and is
        // already unambiguous; only a differing code makes a number collide.
        let mut first = HashMap::new();
        let collides = numbers
            .iter()
            .any(|(_, number, code)| *first.entry(*number).or_insert(code) != code);
        if !collides {
            continue;
        }
        let mut shift = 0;
        let mut previous: Option<(i32, &str)> = None;
        for (residue, number, code) in &numbers {
            if previous.is_some_and(|(n, c)| n == *number && c != code) {
                shift += 1;
            }
            previous = Some((*number, code));
            if shift > 0 {
                renumbered.push((*residue, number.saturating_add(shift)));
            }
        }
    }
    let changed = !renumbered.is_empty();
    for (residue, number) in renumbered {
        let _ = data.topology.residues.set_label_seq_id(residue, number);
    }
    changed
}

#[cfg(test)]
#[path = "chain_labels_tests.rs"]
mod tests;
