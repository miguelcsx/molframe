//! Explicit selection of alternate conformations.
//!
//! Reading never collapses disorder. This module turns a declared analysis
//! policy into an atom selection while leaving the stored rows untouched.

use super::{AtomRef, ResidueRef, Structure};
use crate::contract::{AltlocPolicy, Analysis, AnalysisPolicy, Coverage, Status};
use crate::diagnostic::{Code, Diagnostic};
use crate::selection::AtomSelection;
use crate::symbol::{AltId, SymbolId};
use hashbrown::HashMap;

#[derive(Clone, Copy)]
struct LabelScore {
    label: AltId,
    occupancy: f32,
}

impl Structure {
    /// Selects atom rows under the policy's alternate-conformation rule.
    ///
    /// Blank-labelled atoms are retained under every rule because they belong
    /// to every conformation. The stored structure is never changed.
    #[must_use]
    pub fn resolve_altlocs(&self, policy: &AnalysisPolicy) -> Analysis<AtomSelection> {
        let total = self.atom_count();
        let (selection, missing_label) = match &policy.altloc {
            AltlocPolicy::KeepAll => (AtomSelection::All(total), false),
            AltlocPolicy::ConformerConsistent => (self.conformer_consistent(), false),
            AltlocPolicy::Label(label) => self.named_conformer(label),
            AltlocPolicy::First => (self.per_residue(&FirstLabel), false),
            AltlocPolicy::HighestOccupancyPerResidue => {
                (self.per_residue(&HighestOccupancy), false)
            }
            AltlocPolicy::HighestOccupancyPerAtom => (self.per_atom_occupancy(), false),
        };

        let used = selection.len();
        let mut result = Analysis::complete(
            selection,
            Coverage {
                intended: total,
                used,
                missing: 0,
                ambiguous: total.saturating_sub(used),
            },
            policy,
        );
        if matches!(policy.altloc, AltlocPolicy::KeepAll) && self.has_altlocs() {
            result.status = Status::Ambiguous;
        }
        if policy.altloc.is_hazardous() && self.has_altlocs() {
            result.warnings.push(Diagnostic::new(Code::W6001));
        }
        if missing_label {
            result.warnings.push(
                Diagnostic::new(Code::W4003)
                    .with_context("altloc", named_policy_label(&policy.altloc)),
            );
        }
        result
    }

    fn has_altlocs(&self) -> bool {
        self.data()
            .atoms()
            .any(|atom| atom.alt_id().is_some_and(|label| !label.is_blank()))
    }

    fn named_conformer(&self, label: &str) -> (AtomSelection, bool) {
        let Some(symbol) = self.data().dictionary.get(label) else {
            return (self.select_label(None), true);
        };
        let wanted = AltId::labelled(symbol);
        let present = self
            .data()
            .atoms()
            .any(|atom| atom.alt_id() == Some(wanted));
        let selected = self.select_label(Some(wanted));
        (selected, !present)
    }

    fn select_label(&self, selected: Option<AltId>) -> AtomSelection {
        let positions = self
            .data()
            .atoms()
            .filter(|atom| {
                atom.alt_id()
                    .is_some_and(|label| label.is_blank() || Some(label) == selected)
            })
            .map(|atom| atom.index().get())
            .collect();
        AtomSelection::from_sorted(positions)
    }

    fn conformer_consistent(&self) -> AtomSelection {
        let mut positions = Vec::new();
        for chain in self.data().chains() {
            let chosen = best_label(chain.residues().flat_map(ResidueRef::atoms));
            for atom in chain.residues().flat_map(ResidueRef::atoms) {
                if keeps_label(atom, chosen) {
                    positions.push(atom.index().get());
                }
            }
        }
        AtomSelection::from_sorted(positions)
    }

    fn per_residue(&self, chooser: &impl LabelChooser) -> AtomSelection {
        let mut positions = Vec::new();
        for residue in self.data().residues() {
            let chosen = chooser.choose(residue);
            positions.extend(
                residue
                    .atoms()
                    .filter(|atom| keeps_label(*atom, chosen))
                    .map(|atom| atom.index().get()),
            );
        }
        AtomSelection::from_sorted(positions)
    }

    fn per_atom_occupancy(&self) -> AtomSelection {
        let mut positions = Vec::new();
        for residue in self.data().residues() {
            let mut choices: HashMap<SymbolId, (u32, f32)> = HashMap::new();
            for atom in residue.atoms() {
                let Some(label) = atom.alt_id() else {
                    continue;
                };
                if label.is_blank() {
                    positions.push(atom.index().get());
                    continue;
                }
                let Some(name) = atom.name_symbol() else {
                    continue;
                };
                let occupancy = occupancy(atom);
                choices
                    .entry(name)
                    .and_modify(|choice| {
                        if occupancy > choice.1 {
                            *choice = (atom.index().get(), occupancy);
                        }
                    })
                    .or_insert((atom.index().get(), occupancy));
            }
            positions.extend(choices.into_values().map(|choice| choice.0));
        }
        positions.sort_unstable();
        AtomSelection::from_sorted(positions)
    }
}

trait LabelChooser {
    fn choose(&self, residue: ResidueRef<'_>) -> Option<AltId>;
}

struct FirstLabel;

impl LabelChooser for FirstLabel {
    fn choose(&self, residue: ResidueRef<'_>) -> Option<AltId> {
        residue
            .atoms()
            .filter_map(AtomRef::alt_id)
            .find(|label| !label.is_blank())
    }
}

struct HighestOccupancy;

impl LabelChooser for HighestOccupancy {
    fn choose(&self, residue: ResidueRef<'_>) -> Option<AltId> {
        best_label(residue.atoms())
    }
}

fn best_label<'a>(atoms: impl Iterator<Item = AtomRef<'a>>) -> Option<AltId> {
    let mut scores: Vec<LabelScore> = Vec::new();
    for atom in atoms {
        let Some(label) = atom.alt_id() else {
            continue;
        };
        if label.is_blank() {
            continue;
        }
        let occupancy = occupancy(atom);
        if let Some(score) = scores.iter_mut().find(|score| score.label == label) {
            score.occupancy += occupancy;
        } else {
            scores.push(LabelScore { label, occupancy });
        }
    }
    let mut chosen: Option<LabelScore> = None;
    for score in scores {
        if chosen.is_none_or(|current| score.occupancy > current.occupancy) {
            chosen = Some(score);
        }
    }
    chosen.map(|score| score.label)
}

fn occupancy(atom: AtomRef<'_>) -> f32 {
    match atom.occupancy() {
        Some(value) => value,
        None => 1.0,
    }
}

fn keeps_label(atom: AtomRef<'_>, selected: Option<AltId>) -> bool {
    atom.alt_id()
        .is_some_and(|label| label.is_blank() || Some(label) == selected)
}

fn named_policy_label(policy: &AltlocPolicy) -> Box<str> {
    match policy {
        AltlocPolicy::Label(label) => label.clone(),
        _ => "".into(),
    }
}

#[cfg(test)]
#[path = "disorder_tests.rs"]
mod tests;
