//! The raw component tables a structure file carries, indexed without lowering.
//!
//! A [`Source`] keeps the `chem_comp`, `chem_comp_atom` and `chem_comp_bond`
//! categories of one data block as read, plus two sorted run tables that say
//! which rows belong to which component. Building the tables is one pass over
//! the `comp_id` columns, `O(rows + k log k)` for `k` runs, with one allocation
//! per table and no per-row allocation: ids are shared with the category's own
//! text. Checking a source reads every row through the borrowed readers of
//! `views` and copies nothing; only a source that fails the check is lowered
//! whole, to report exactly what the eager lowering would have reported.

use super::lowering::{Rows, component_from};
use super::views::{AtomView, BondView, atom_view, bond_view, coordinate_rows, topology_is_valid};
use crate::Component;
use molframe_cif::{Category, CifValue};
use molframe_core::diagnostic::{Code, Diagnostic};
use std::collections::BTreeSet;
use std::sync::{Arc, OnceLock};

/// Consecutive rows of one category that share a component id.
#[derive(Clone, Debug)]
pub(super) struct Run {
    id: Arc<str>,
    start: usize,
    end: usize,
}

/// A component the source defines: its `chem_comp` row and, once asked for,
/// its lowered form.
#[derive(Clone, Debug)]
pub(super) struct Entry {
    pub(super) id: Arc<str>,
    row: usize,
    resolved: OnceLock<Result<Arc<Component>, Diagnostic>>,
}

impl Entry {
    pub(super) fn component(&self, source: &Source) -> Result<Arc<Component>, Diagnostic> {
        self.resolved
            .get_or_init(|| source.lower(self.row, &self.id).map(Arc::new))
            .clone()
    }
}

/// One data block's component tables and the index over them.
#[derive(Clone, Debug)]
pub(super) struct Source {
    metadata: Category,
    atoms: Category,
    bonds: Option<Category>,
    atom_runs: Vec<Run>,
    bond_runs: Vec<Run>,
    /// Defined components, sorted by id.
    pub(super) entries: Vec<Entry>,
    /// Whether the atom rows state formal charges.
    pub(super) charged: bool,
}

impl Source {
    /// Indexes the tables; `None` when the block defines no component, which
    /// is the usual case for a deposited file that only lists them.
    pub(super) fn new(
        metadata: Category,
        atoms: Category,
        bonds: Option<Category>,
    ) -> Option<Self> {
        let atom_runs = runs(&atoms);
        let bond_runs = bonds.as_ref().map(runs).unwrap_or_default();
        let mut entries = Vec::new();
        for row in 0..metadata.row_count() {
            let Some(id) = shared_id(&metadata, "id", row) else {
                continue;
            };
            if !runs_of(&atom_runs, &id).is_empty() {
                entries.push(Entry {
                    id,
                    row,
                    resolved: OnceLock::new(),
                });
            }
        }
        if entries.is_empty() {
            return None;
        }
        entries.sort_by(|left, right| left.id.cmp(&right.id).then(left.row.cmp(&right.row)));
        let charged = atoms.column("charge").is_some();
        Some(Self {
            metadata,
            atoms,
            bonds,
            atom_runs,
            bond_runs,
            entries,
            charged,
        })
    }

    /// The entry for `id`, by binary search.
    pub(super) fn entry(&self, id: &str) -> Option<&Entry> {
        let at = self.entries.partition_point(|entry| entry.id.as_ref() < id);
        self.entries.get(at).filter(|entry| entry.id.as_ref() == id)
    }

    fn lower(&self, row: usize, id: &str) -> Result<Component, Diagnostic> {
        let atom_rows: Vec<usize> = rows_of(&self.atom_runs, id).collect();
        let bond_rows: Vec<usize> = rows_of(&self.bond_runs, id).collect();
        let bonds = self
            .bonds
            .as_ref()
            .filter(|_| !bond_rows.is_empty())
            .map(|category| Rows {
                category,
                rows: &bond_rows,
            });
        component_from(
            &self.metadata,
            row,
            Some(Rows {
                category: &self.atoms,
                rows: &atom_rows,
            }),
            bonds,
            true,
        )
        .map_err(|found| {
            found
                .into_iter()
                .next()
                .unwrap_or_else(|| Diagnostic::new(Code::W3206).in_category(self.atoms.name()))
        })
    }

    /// Whether every defined component lowers, decided without copying a row.
    ///
    /// Reads each row of each component once through the same readers the
    /// lowering uses, so it accepts exactly what the lowering accepts. Time is
    /// `O(rows)` plus the sort of each component's atoms and bonds; the two
    /// scratch tables grow to the largest component and are reused.
    pub(super) fn is_valid(&self) -> bool {
        let mut atoms: Vec<AtomView<'_>> = Vec::new();
        let mut bonds: Vec<BondView<'_>> = Vec::new();
        self.entries
            .iter()
            .all(|entry| self.entry_is_valid(entry, &mut atoms, &mut bonds))
    }

    fn entry_is_valid<'a>(
        &'a self,
        entry: &Entry,
        atoms: &mut Vec<AtomView<'a>>,
        bonds: &mut Vec<BondView<'a>>,
    ) -> bool {
        let id = entry.id.as_ref();
        let named = |item: &str| self.metadata.identifier(item, entry.row).is_some();
        if !named("name") || self.metadata.text("type", entry.row).is_none() {
            return false;
        }
        atoms.clear();
        bonds.clear();
        for row in rows_of(&self.atom_runs, id) {
            match atom_view(&self.atoms, row, true) {
                Ok(view) => atoms.push(view),
                Err(_) => return false,
            }
        }
        if let Some(category) = &self.bonds {
            for row in rows_of(&self.bond_runs, id) {
                match bond_view(category, row) {
                    Ok(view) => bonds.push(view),
                    Err(_) => return false,
                }
            }
        }
        if !topology_is_valid(atoms, bonds) {
            return false;
        }
        [true, false].into_iter().all(|ideal| {
            coordinate_rows(&self.atoms, rows_of(&self.atom_runs, id), ideal, drop).is_ok()
        })
    }

    /// Every finding the lowering of this source reports, in file order.
    pub(super) fn findings(&self) -> Vec<Diagnostic> {
        let mut findings = Vec::new();
        for row in 0..self.metadata.row_count() {
            let Some(id) = shared_id(&self.metadata, "id", row) else {
                continue;
            };
            if runs_of(&self.atom_runs, &id).is_empty() {
                continue;
            }
            let atom_rows: Vec<usize> = rows_of(&self.atom_runs, &id).collect();
            let bond_rows: Vec<usize> = rows_of(&self.bond_runs, &id).collect();
            findings.extend(self.report(row, &atom_rows, &bond_rows));
        }
        findings
    }

    /// The first component id this source defines that `seen` already holds,
    /// in file order; every id it defines is added to `seen`.
    pub(super) fn first_duplicate(&self, seen: &mut BTreeSet<Arc<str>>) -> Option<Arc<str>> {
        let mut duplicate = None;
        for row in 0..self.metadata.row_count() {
            let Some(id) = shared_id(&self.metadata, "id", row) else {
                continue;
            };
            if !runs_of(&self.atom_runs, &id).is_empty()
                && !seen.insert(Arc::clone(&id))
                && duplicate.is_none()
            {
                duplicate = Some(id);
            }
        }
        duplicate
    }

    fn report(&self, row: usize, atom_rows: &[usize], bond_rows: &[usize]) -> Vec<Diagnostic> {
        let bonds = self.bonds.as_ref().filter(|_| !bond_rows.is_empty());
        component_from(
            &self.metadata,
            row,
            Some(Rows {
                category: &self.atoms,
                rows: atom_rows,
            }),
            bonds.map(|category| Rows {
                category,
                rows: bond_rows,
            }),
            true,
        )
        .err()
        .unwrap_or_default()
    }
}

/// The id in one row as shared text, copying only a number written as one.
fn shared_id(category: &Category, item: &str, row: usize) -> Option<Arc<str>> {
    match category.value(item, row)? {
        CifValue::Text(text) => Some(Arc::clone(text)),
        other => other.as_identifier().map(|id| Arc::from(id.as_ref())),
    }
}

/// Maximal runs of equal `comp_id`, sorted by id then position.
fn runs(category: &Category) -> Vec<Run> {
    let mut runs: Vec<Run> = Vec::new();
    for row in 0..category.row_count() {
        let Some(id) = shared_id(category, "comp_id", row) else {
            continue;
        };
        match runs.last_mut() {
            Some(last) if last.end == row && last.id == id => last.end = row + 1,
            _ => runs.push(Run {
                id,
                start: row,
                end: row + 1,
            }),
        }
    }
    runs.sort_by(|left, right| left.id.cmp(&right.id).then(left.start.cmp(&right.start)));
    runs
}

/// The runs of one id, found by binary search.
fn runs_of<'a>(runs: &'a [Run], id: &str) -> &'a [Run] {
    let start = runs.partition_point(|run| run.id.as_ref() < id);
    let len = runs[start..]
        .iter()
        .take_while(|run| run.id.as_ref() == id)
        .count();
    &runs[start..start + len]
}

/// Every row of one id, in file order, without collecting them.
fn rows_of<'a>(runs: &'a [Run], id: &str) -> impl Iterator<Item = usize> + Clone + use<'a> {
    runs_of(runs, id).iter().flat_map(|run| run.start..run.end)
}
