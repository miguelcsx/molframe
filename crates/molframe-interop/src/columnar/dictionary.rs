//! Shared Arrow dictionaries for categorical atom-table text columns.
//!
//! A structure has a handful of distinct chains, residue names, atom names and
//! elements, however many atoms it holds. Each categorical column is therefore
//! exported as an Arrow `Dictionary(Int32, Utf8)` whose values array is built
//! once per table and shared (reference-counted, never copied) by every record
//! batch. Per atom the export holds one `i32` key and never a string.
//!
//! Keys are assigned in ascending raw-identifier order over the identifiers
//! that actually occur, so the dictionary has no unused categories and its
//! numbering is deterministic.

use crate::numeric::usize_to_u32;
use arrow::array::{ArrayRef, DictionaryArray, Int32Array, StringArray};
use arrow::datatypes::{DataType, Int32Type};
use arrow::error::Result;
use molframe_core::topology::{ChainTable, ResidueTable};
use molframe_core::{AltId, ChainIndex, Element, ResidueIndex, Structure, SymbolId};
use std::sync::Arc;

/// Key stored for an atom whose value is absent; the exported key is null.
const ABSENT: i32 = -1;

/// The Arrow type of every categorical text column.
pub(crate) fn dictionary_type() -> DataType {
    DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8))
}

/// Maps raw identifiers to compact dictionary keys.
#[derive(Debug)]
pub(crate) struct SymbolDictionary {
    keys: Vec<i32>,
    values: ArrayRef,
}

impl SymbolDictionary {
    /// Compact key of `raw`, or `None` when the identifier does not occur.
    pub(crate) fn key(&self, raw: u32) -> Option<i32> {
        self.keys
            .get(raw as usize)
            .copied()
            .filter(|key| *key != ABSENT)
    }

    /// Number of distinct values.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.values.len()
    }

    /// Wraps per-row keys with the shared values array. `O(rows)`.
    pub(crate) fn array(&self, keys: impl IntoIterator<Item = Option<i32>>) -> Result<ArrayRef> {
        let keys: Int32Array = keys.into_iter().collect();
        Ok(Arc::new(DictionaryArray::<Int32Type>::try_new(
            keys,
            self.values.clone(),
        )?))
    }
}

/// Collects the identifiers that occur, then fixes the dictionary.
#[derive(Default)]
struct Occurrences {
    used: Vec<bool>,
}

impl Occurrences {
    fn observe(&mut self, raw: u32) {
        let raw = raw as usize;
        if self.used.len() <= raw {
            self.used.resize(raw + 1, false);
        }
        self.used[raw] = true;
    }

    fn finish<'a>(self, mut text: impl FnMut(u32) -> Option<&'a str>) -> SymbolDictionary {
        let mut keys = vec![ABSENT; self.used.len()];
        let mut values = Vec::new();
        for (raw, used) in self.used.iter().enumerate() {
            if *used {
                // Compact keys are bounded by the number of distinct strings.
                keys[raw] = match i32::try_from(values.len()) {
                    Ok(key) => key,
                    Err(_) => i32::MAX,
                };
                values.push(match text(usize_to_u32(raw)) {
                    Some(text) => text,
                    None => "",
                });
            }
        }
        SymbolDictionary {
            keys,
            values: Arc::new(StringArray::from(values)),
        }
    }
}

/// Dictionaries for every categorical atom column of one structure.
///
/// Built once, lazily, in `O(atoms)` time and `O(distinct values + residues)`
/// space, then shared by every batch.
#[derive(Debug)]
pub(crate) struct AtomDictionaries {
    pub(crate) atom_name: SymbolDictionary,
    pub(crate) altloc: SymbolDictionary,
    pub(crate) element: SymbolDictionary,
    pub(crate) chain: SymbolDictionary,
    pub(crate) label_asym: SymbolDictionary,
    pub(crate) resname: SymbolDictionary,
    /// Per-residue key into `chain`.
    pub(crate) residue_chain: Vec<i32>,
    /// Per-residue key into `label_asym`.
    pub(crate) residue_label_asym: Vec<i32>,
    /// Per-residue key into `resname`.
    pub(crate) residue_resname: Vec<i32>,
}

impl AtomDictionaries {
    pub(crate) fn build(structure: &Structure) -> Self {
        let data = structure.data();
        let text = |raw: u32| data.dictionary.resolve(SymbolId::from_raw(raw));

        let mut names = Occurrences::default();
        let mut altlocs = Occurrences::default();
        let mut elements = Occurrences::default();
        for chunk in data.chunks.iter() {
            match chunk.atom_names_plain() {
                Some(plain) => plain.iter().for_each(|symbol| names.observe(symbol.get())),
                None => (0..chunk.len())
                    .filter_map(|local| chunk.atom_name(local))
                    .for_each(|symbol| names.observe(symbol.get())),
            }
            for local in 0..chunk.len() {
                if let Some(alt) = chunk.alt_id(local)
                    && !alt.is_blank()
                {
                    altlocs.observe(alt.get());
                }
                if let Some(element) = chunk.element(local) {
                    elements.observe(u32::from(element.atomic_number()));
                }
            }
        }

        let (chain, label_asym, residue_chain, residue_label_asym) =
            chain_dictionaries(&data.topology.chains, &data.topology.residues, structure);
        let (resname, residue_resname) = resname_dictionary(&data.topology.residues, structure);

        Self {
            atom_name: names.finish(text),
            altloc: altlocs.finish(|raw| {
                AltId::from_raw(raw)
                    .symbol()
                    .and_then(|symbol| data.dictionary.resolve(symbol))
            }),
            element: elements.finish(|raw| match u8::try_from(raw) {
                Ok(number) => Some(Element::from_atomic_number(number).symbol()),
                Err(_) => None,
            }),
            chain,
            label_asym,
            resname,
            residue_chain,
            residue_label_asym,
            residue_resname,
        }
    }
}

type ChainDictionaries = (SymbolDictionary, SymbolDictionary, Vec<i32>, Vec<i32>);

/// Author chain identifier (falling back to the label one) and label chain
/// identifier, each with a per-residue key vector. `O(chains + residues)`.
fn chain_dictionaries(
    chains: &ChainTable,
    residues: &ResidueTable,
    structure: &Structure,
) -> ChainDictionaries {
    let data = structure.data();
    let text = |raw: u32| data.dictionary.resolve(SymbolId::from_raw(raw));
    let identifiers = (0..chains.len())
        .map(|position| {
            let chain = ChainIndex::new(usize_to_u32(position));
            let label = chains.label_asym_id(chain);
            (chain, chains.auth_asym_id(chain).or(label), label)
        })
        .collect::<Vec<_>>();
    let mut author = Occurrences::default();
    let mut label = Occurrences::default();
    for (_, author_id, label_id) in &identifiers {
        if let Some(symbol) = author_id {
            author.observe(symbol.get());
        }
        if let Some(symbol) = label_id {
            label.observe(symbol.get());
        }
    }
    let author = author.finish(text);
    let label = label.finish(text);
    let mut residue_author = vec![ABSENT; residues.len()];
    let mut residue_label = vec![ABSENT; residues.len()];
    for (chain, author_id, label_id) in &identifiers {
        let Some(range) = chains.residues(*chain) else {
            continue;
        };
        let author_key = author_id.and_then(|symbol| author.key(symbol.get()));
        let label_key = label_id.and_then(|symbol| label.key(symbol.get()));
        for residue in range {
            if let Some(slot) = residue_author.get_mut(residue as usize) {
                *slot = key_or_absent(author_key);
            }
            if let Some(slot) = residue_label.get_mut(residue as usize) {
                *slot = key_or_absent(label_key);
            }
        }
    }
    (author, label, residue_author, residue_label)
}

/// Residue-name dictionary and a per-residue key vector. `O(residues)`.
fn resname_dictionary(
    residues: &ResidueTable,
    structure: &Structure,
) -> (SymbolDictionary, Vec<i32>) {
    let data = structure.data();
    let symbols = (0..residues.len())
        .map(|position| {
            let residue = ResidueIndex::new(usize_to_u32(position));
            residues
                .label_comp_id(residue)
                .or_else(|| residues.auth_comp_id(residue))
        })
        .collect::<Vec<_>>();
    let mut occurrences = Occurrences::default();
    symbols
        .iter()
        .flatten()
        .for_each(|symbol| occurrences.observe(symbol.get()));
    let dictionary = occurrences.finish(|raw| data.dictionary.resolve(SymbolId::from_raw(raw)));
    let keys = symbols
        .iter()
        .map(|symbol| key_or_absent(symbol.and_then(|symbol| dictionary.key(symbol.get()))))
        .collect();
    (dictionary, keys)
}

fn key_or_absent(key: Option<i32>) -> i32 {
    match key {
        Some(key) => key,
        None => ABSENT,
    }
}

/// `Some(key)` unless the key marks an absent value.
pub(crate) fn present(key: i32) -> Option<i32> {
    (key != ABSENT).then_some(key)
}

#[cfg(test)]
#[path = "dictionary_tests.rs"]
mod tests;
