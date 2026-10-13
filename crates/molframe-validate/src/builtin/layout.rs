//! Residue layout shared by the built-in checks.
//!
//! One sequential pass records, for every residue, which chain it is in and
//! its parent residue name (protonation variants folded onto the standard
//! residue). The checks then run over residue ranges in parallel and ask this
//! layout for neighbours and atom positions, so none of them walks the chain
//! structure again.

use molframe_analysis::parent_residue;
use molframe_core::index::ResidueIndex;
use molframe_core::structure::Structure;

/// The twenty standard amino acids.
const AMINO_ACIDS: [&str; 20] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Longest C-N distance still treated as a peptide link, in ångström.
const PEPTIDE_LINK: f32 = 2.0;
/// Longest O3'-P distance still treated as a phosphodiester link.
const PHOSPHATE_LINK: f32 = 2.2;

/// Neighbour and position lookups over one structure.
#[derive(Debug)]
pub struct Layout<'a> {
    structure: &'a Structure,
    positions: &'a [[f32; 3]],
    chain: Vec<u32>,
    name: Vec<&'a str>,
}

impl<'a> Layout<'a> {
    /// Records chain membership and parent names for every residue.
    ///
    /// Runs in `O(residues)` time.
    #[must_use]
    pub fn new(structure: &'a Structure) -> Self {
        let count = structure.residue_count();
        let mut chain = vec![u32::MAX; count];
        let mut name = vec![""; count];
        for (ordinal, chain_ref) in structure.data().chains().enumerate() {
            let label = match u32::try_from(ordinal) {
                Ok(label) => label,
                Err(_) => u32::MAX,
            };
            for residue in chain_ref.residues() {
                let at = residue.index().as_usize();
                if let (Some(slot), Some(text)) = (chain.get_mut(at), name.get_mut(at)) {
                    *slot = label;
                    *text = match residue.name() {
                        Some(raw) => parent_residue(raw),
                        None => "",
                    };
                }
            }
        }
        Self {
            structure,
            positions: structure.positions(),
            chain,
            name,
        }
    }

    /// Number of residues.
    #[must_use]
    pub fn len(&self) -> usize {
        self.name.len()
    }

    /// Whether the structure has no residues.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.name.is_empty()
    }

    /// The parent residue name at this index, empty when unknown.
    #[must_use]
    pub fn name(&self, residue: usize) -> &'a str {
        match self.name.get(residue) {
            Some(name) => name,
            None => "",
        }
    }

    /// Whether this is one of the twenty standard amino acids.
    #[must_use]
    pub fn is_amino_acid(&self, residue: usize) -> bool {
        AMINO_ACIDS.contains(&self.name(residue))
    }

    /// Whether the residue carries the ring atoms of a nucleotide sugar.
    #[must_use]
    pub fn is_nucleotide(&self, residue: usize) -> bool {
        ["C1'", "C3'", "C4'", "O4'"]
            .iter()
            .all(|atom| self.position(residue, atom).is_some())
            && !self.is_amino_acid(residue)
    }

    /// The atom index and position of a named atom, first alternate location.
    #[must_use]
    pub fn atom(&self, residue: usize, atom: &str) -> Option<(u32, [f32; 3])> {
        let ordinal = u32::try_from(residue).ok()?;
        let handle = self.structure.data().residue(ResidueIndex::new(ordinal))?;
        let index = handle.atom(atom)?.index();
        let position = self.positions.get(index.as_usize()).copied()?;
        Some((index.get(), position))
    }

    /// The position of a named atom in a residue, first alternate location.
    #[must_use]
    pub fn position(&self, residue: usize, atom: &str) -> Option<[f32; 3]> {
        self.atom(residue, atom).map(|(_, position)| position)
    }

    fn same_chain(&self, previous: usize, next: usize) -> bool {
        match (self.chain.get(previous), self.chain.get(next)) {
            (Some(left), Some(right)) => left == right && *left != u32::MAX,
            _ => false,
        }
    }

    /// Whether `residue` is joined to the one before it by a peptide bond.
    #[must_use]
    pub fn peptide_linked(&self, residue: usize) -> bool {
        self.linked(residue, ("C", "N"), PEPTIDE_LINK)
    }

    /// Whether `residue` is joined to the one before it by a phosphodiester.
    #[must_use]
    pub fn phosphate_linked(&self, residue: usize) -> bool {
        self.linked(residue, ("O3'", "P"), PHOSPHATE_LINK)
    }

    fn linked(&self, residue: usize, atoms: (&str, &str), limit: f32) -> bool {
        let Some(previous) = residue.checked_sub(1) else {
            return false;
        };
        if !self.same_chain(previous, residue) {
            return false;
        }
        match (
            self.position(previous, atoms.0),
            self.position(residue, atoms.1),
        ) {
            (Some(left), Some(right)) => molframe_geom::distance(left, right) <= f64::from(limit),
            _ => false,
        }
    }
}

#[cfg(test)]
#[path = "layout_tests.rs"]
mod tests;
