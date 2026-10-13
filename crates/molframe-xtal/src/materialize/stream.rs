//! Memory-bounded expansion of an assembly into consecutive structures.

use super::{LabelLedger, Materializer, capacity, single};
use crate::AssemblyView;
use molframe_core::{Diagnostic, Structure};

impl AssemblyView {
    /// Expands the assembly in bounded pieces instead of all at once.
    ///
    /// Each yielded structure holds whole placements (every chain copy that
    /// shares one transform), so bonds between chains of one placement survive
    /// and no placement is split. A piece closes once it holds at least
    /// `max_atoms` atoms; one placement larger than the bound is yielded alone.
    /// Chain labels are collision-free across the whole sequence and equal to
    /// those of [`AssemblyView::materialize`], atom-site identifiers continue
    /// from piece to piece, and every atom carries its assembly-wide instance
    /// identifier. Only one piece is alive at a time, so peak memory follows
    /// `max_atoms`, not the size of the assembly.
    #[must_use]
    pub fn chunks(&self, max_atoms: usize) -> AssemblyChunks<'_> {
        AssemblyChunks {
            view: self,
            max_atoms: max_atoms.max(1),
            cursor: 0,
            atom_base: 0,
            labels: Some(LabelLedger::default()),
        }
    }
}

/// Iterator over the pieces of an expanded assembly; see [`AssemblyView::chunks`].
#[derive(Debug)]
pub struct AssemblyChunks<'a> {
    view: &'a AssemblyView,
    max_atoms: usize,
    cursor: usize,
    atom_base: u32,
    labels: Option<LabelLedger>,
}

impl AssemblyChunks<'_> {
    /// The end of the next piece: whole transform runs up to the atom bound.
    fn next_end(&self) -> usize {
        let instances = &self.view.instances;
        let mut end = self.cursor;
        let mut atoms = 0usize;
        while end < instances.len() {
            let transform = instances[end].transform;
            let run_start = end;
            while end < instances.len() && instances[end].transform == transform {
                atoms = atoms.saturating_add(instances[end].atoms.len());
                end += 1;
            }
            if atoms >= self.max_atoms && run_start < end {
                break;
            }
        }
        end
    }

    fn expand(&mut self, end: usize) -> Result<Structure, Vec<Diagnostic>> {
        let piece = self.view.slice(self.cursor..end);
        let labels = self.labels.take().unwrap_or_default();
        let (structure, labels) = Materializer::new(&piece, labels, self.atom_base)
            .map_err(single)?
            .build()?;
        let added = structure.atom_count();
        self.atom_base = self
            .atom_base
            .checked_add(added)
            .ok_or_else(|| single(capacity("atoms")))?;
        self.labels = Some(labels);
        Ok(structure)
    }
}

impl Iterator for AssemblyChunks<'_> {
    type Item = Result<Structure, Vec<Diagnostic>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.view.instances.len() {
            return None;
        }
        let end = self.next_end();
        let result = self.expand(end);
        // A failed piece ends the sequence: the ledger is gone and any later
        // piece would carry inconsistent labels.
        self.cursor = if result.is_ok() {
            end
        } else {
            self.view.instances.len()
        };
        Some(result)
    }
}
