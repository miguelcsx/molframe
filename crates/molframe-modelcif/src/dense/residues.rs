//! Residue identities a dense confidence matrix is indexed by.
//!
//! Chain and component names are interned once, so each residue costs sixteen
//! bytes regardless of how long its labels are.

use molframe_core::numeric::usize_to_u32;

const NO_COMPONENT: u32 = u32::MAX;

/// One residue named by its label chain, label sequence number and component.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ResidueRef<'a> {
    /// Label chain identifier.
    pub chain: &'a str,
    /// Label sequence number.
    pub sequence: i64,
    /// Chemical component identifier, when the source stated one.
    pub component: Option<&'a str>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Entry {
    chain: u32,
    component: u32,
    sequence: i64,
}

/// Ordered residues of a dense matrix; position `i` labels row and column `i`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ResidueIndex {
    chains: Vec<Box<str>>,
    components: Vec<Box<str>>,
    entries: Vec<Entry>,
}

impl ResidueIndex {
    /// Creates an empty index.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            chains: Vec::new(),
            components: Vec::new(),
            entries: Vec::new(),
        }
    }

    /// Appends a residue and returns its position.
    pub fn push(&mut self, chain: &str, sequence: i64, component: Option<&str>) -> usize {
        let chain = intern(&mut self.chains, chain);
        let component = component.map_or(NO_COMPONENT, |name| intern(&mut self.components, name));
        self.entries.push(Entry {
            chain,
            component,
            sequence,
        });
        self.entries.len() - 1
    }

    /// Number of residues.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the index has no residues.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Residue at `position`.
    #[must_use]
    pub fn get(&self, position: usize) -> Option<ResidueRef<'_>> {
        let entry = self.entries.get(position)?;
        Some(ResidueRef {
            chain: self.chains.get(entry.chain as usize)?,
            sequence: entry.sequence,
            component: self
                .components
                .get(entry.component as usize)
                .map(AsRef::as_ref),
        })
    }

    /// Position of the residue with this chain and sequence number.
    #[must_use]
    pub fn position(&self, chain: &str, sequence: i64) -> Option<usize> {
        let chain = self.chains.iter().position(|name| name.as_ref() == chain)?;
        self.entries
            .iter()
            .position(|entry| entry.chain as usize == chain && entry.sequence == sequence)
    }

    /// Iterates residues in matrix order.
    pub fn iter(&self) -> impl Iterator<Item = ResidueRef<'_>> {
        (0..self.entries.len()).filter_map(|position| self.get(position))
    }

    /// Whether any residue states a component.
    #[must_use]
    pub fn has_components(&self) -> bool {
        self.entries
            .iter()
            .any(|entry| entry.component != NO_COMPONENT)
    }
}

fn intern(names: &mut Vec<Box<str>>, name: &str) -> u32 {
    // Residues arrive chain by chain, so the newest names are the likeliest hits.
    if let Some(position) = names
        .iter()
        .rposition(|candidate| candidate.as_ref() == name)
    {
        return usize_to_u32(position);
    }
    names.push(name.into());
    usize_to_u32(names.len() - 1)
}

#[cfg(test)]
#[path = "residues_tests.rs"]
mod tests;
