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

/// Accumulated occupancy for one alternate-conformation label.
///
/// `order` records when the label was first encountered so occupancy ties
/// preserve the original deterministic first-label behavior.
#[derive(Clone, Copy)]
struct LabelScore {
    occupancy: f32,
    order: usize,
}

/// The label selected for one residue or chain and whether disorder was seen.
#[derive(Clone, Copy)]
struct LabelChoice {
    selected: Option<AltId>,
    has_altlocs: bool,
}

/// Internal result of resolving one alternate-conformation policy.
///
/// Besides the selected atoms, resolution carries metadata already discovered
/// while traversing atoms so callers do not need additional full-structure
/// scans.
struct AltlocResolution {
    selection: AtomSelection,
    missing_label: bool,
    has_altlocs: bool,
}

/// Best retained row for one atom name under per-atom occupancy resolution.
#[derive(Clone, Copy)]
struct AtomChoice {
    atom: u32,
    occupancy: f32,
}

impl Structure {
    /// Selects atom rows under the policy's alternate-conformation rule.
    ///
    /// Blank-labelled atoms are retained under every rule because they belong
    /// to every conformation. The stored structure is never changed.
    ///
    /// # Panics
    ///
    /// This function does not intentionally panic when an internal selector
    /// reports an impossible cardinality; such a result is returned as an
    /// indeterminate analysis instead.
    #[must_use]
    pub fn resolve_altlocs(&self, policy: &AnalysisPolicy) -> Analysis<AtomSelection> {
        let total = self.atom_count();

        let AltlocResolution {
            selection,
            missing_label,
            has_altlocs,
        } = match &policy.altloc {
            AltlocPolicy::KeepAll => AltlocResolution {
                selection: AtomSelection::All(total),
                missing_label: false,
                has_altlocs: self.has_altlocs(),
            },
            AltlocPolicy::ConformerConsistent => self.conformer_consistent(),
            AltlocPolicy::Label(label) => self.named_conformer(label),
            AltlocPolicy::First => self.per_residue(&FirstLabel),
            AltlocPolicy::HighestOccupancyPerResidue => self.per_residue(&HighestOccupancy),
            AltlocPolicy::HighestOccupancyPerAtom => self.per_atom_occupancy(),
        };

        let selected_rows = selection.len();
        let used = match u32::try_from(selected_rows) {
            Ok(used) if used <= total => used,
            _ => {
                return Analysis::indeterminate(
                    selection,
                    Coverage {
                        intended: total,
                        used: 0,
                        missing: total,
                        ambiguous: 0,
                    },
                    policy,
                )
                .with_warning(
                    Diagnostic::new(Code::E3001)
                        .with_context("selection rows", selected_rows.to_string())
                        .with_context("atoms", total.to_string()),
                );
            }
        };

        let mut result = Analysis::complete(
            selection,
            Coverage {
                intended: total,
                used,
                missing: 0,
                ambiguous: total - used,
            },
            policy,
        );

        if matches!(&policy.altloc, AltlocPolicy::KeepAll) && has_altlocs {
            result.status = Status::Ambiguous;
        }

        if policy.altloc.is_hazardous() && has_altlocs {
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

    /// Returns whether the structure contains any non-blank alternate label.
    ///
    /// This scan is needed only for policies such as `KeepAll` that otherwise
    /// do not inspect individual atom rows.
    fn has_altlocs(&self) -> bool {
        self.data()
            .atoms()
            .any(|atom| atom.alt_id().is_some_and(|label| !label.is_blank()))
    }

    /// Resolves a specifically named alternate conformation.
    ///
    /// Blank atoms are always retained. A label absent from either the symbol
    /// dictionary or the atom rows yields the blank-only selection and marks
    /// the requested label as missing.
    fn named_conformer(&self, label: &str) -> AltlocResolution {
        let wanted = self.data().dictionary.get(label).and_then(AltId::labelled);

        let requested_valid = wanted.is_some();
        let (selection, present, has_altlocs) = self.select_label(wanted);

        AltlocResolution {
            selection,
            missing_label: !requested_valid || !present,
            has_altlocs,
        }
    }

    /// Selects blank rows together with one optional alternate label.
    ///
    /// Returns the selection, whether the requested non-blank label was
    /// encountered, and whether any alternate conformations exist.
    fn select_label(&self, selected: Option<AltId>) -> (AtomSelection, bool, bool) {
        let mut positions = Vec::new();
        let mut selected_present = false;
        let mut has_altlocs = false;

        for atom in self.data().atoms() {
            let Some(label) = atom.alt_id() else {
                continue;
            };

            if !label.is_blank() {
                has_altlocs = true;
            }

            if Some(label) == selected {
                selected_present = true;
            }

            if label.is_blank() || Some(label) == selected {
                positions.push(atom.index().get());
            }
        }

        (
            AtomSelection::from_sorted(positions),
            selected_present,
            has_altlocs,
        )
    }

    /// Chooses one occupancy-maximizing label consistently for each chain.
    ///
    /// A reusable score table avoids allocating fresh aggregation storage for
    /// every chain while preserving deterministic first-seen tie breaking.
    fn conformer_consistent(&self) -> AltlocResolution {
        let mut positions = Vec::new();
        let mut scores = HashMap::new();
        let mut has_altlocs = false;

        for chain in self.data().chains() {
            let choice = best_label(chain.residues().flat_map(ResidueRef::atoms), &mut scores);

            has_altlocs |= choice.has_altlocs;

            for atom in chain.residues().flat_map(ResidueRef::atoms) {
                if keeps_label(atom, choice.selected) {
                    positions.push(atom.index().get());
                }
            }
        }

        AltlocResolution {
            selection: AtomSelection::from_sorted(positions),
            missing_label: false,
            has_altlocs,
        }
    }

    /// Resolves alternate conformations independently for every residue.
    ///
    /// The chooser receives reusable occupancy scratch storage so repeated
    /// residue resolution retains hash-table capacity instead of reallocating.
    fn per_residue(&self, chooser: &impl LabelChooser) -> AltlocResolution {
        let mut positions = Vec::new();
        let mut scores = HashMap::new();
        let mut has_altlocs = false;

        for residue in self.data().residues() {
            let choice = chooser.choose(residue, &mut scores);
            has_altlocs |= choice.has_altlocs;

            positions.extend(
                residue
                    .atoms()
                    .filter(|atom| keeps_label(*atom, choice.selected))
                    .map(|atom| atom.index().get()),
            );
        }

        AltlocResolution {
            selection: AtomSelection::from_sorted(positions),
            missing_label: false,
            has_altlocs,
        }
    }

    /// Selects the highest-occupancy alternate row independently per atom name.
    ///
    /// Choices are reused across residues and selected atoms are emitted during
    /// a second ordered residue traversal, avoiding a global sort of results.
    fn per_atom_occupancy(&self) -> AltlocResolution {
        let mut positions = Vec::new();
        let mut choices: HashMap<SymbolId, AtomChoice> = HashMap::new();
        let mut has_altlocs = false;

        for residue in self.data().residues() {
            choices.clear();

            for atom in residue.atoms() {
                let Some(label) = atom.alt_id() else {
                    continue;
                };

                if label.is_blank() {
                    continue;
                }

                has_altlocs = true;

                let Some(name) = atom.name_symbol() else {
                    continue;
                };

                record_atom_choice(&mut choices, name, atom);
            }

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

                let Some(choice) = choices.get(&name) else {
                    continue;
                };

                if choice.atom == atom.index().get() {
                    positions.push(choice.atom);
                }
            }
        }

        AltlocResolution {
            selection: AtomSelection::from_sorted(positions),
            missing_label: false,
            has_altlocs,
        }
    }
}

/// Chooses one alternate label for a residue.
///
/// Implementations may use the supplied score map as reusable scratch storage;
/// callers clear or overwrite it through the chooser implementation.
trait LabelChooser {
    /// Returns the selected label together with whether the residue is disordered.
    fn choose(
        &self,
        residue: ResidueRef<'_>,
        scores: &mut HashMap<AltId, LabelScore>,
    ) -> LabelChoice;
}

/// Chooses the first non-blank alternate label in atom order.
struct FirstLabel;

impl LabelChooser for FirstLabel {
    /// Returns the first non-blank alternate label present in the residue.
    fn choose(
        &self,
        residue: ResidueRef<'_>,
        _scores: &mut HashMap<AltId, LabelScore>,
    ) -> LabelChoice {
        let selected = residue
            .atoms()
            .filter_map(AtomRef::alt_id)
            .find(|label| !label.is_blank());

        LabelChoice {
            selected,
            has_altlocs: selected.is_some(),
        }
    }
}

/// Chooses the alternate label with the greatest accumulated occupancy.
struct HighestOccupancy;

impl LabelChooser for HighestOccupancy {
    /// Aggregates occupancy by label and preserves first-seen ordering on ties.
    fn choose(
        &self,
        residue: ResidueRef<'_>,
        scores: &mut HashMap<AltId, LabelScore>,
    ) -> LabelChoice {
        best_label(residue.atoms(), scores)
    }
}

/// Finds the occupancy-maximizing non-blank alternate label.
///
/// The reusable hash table gives expected constant-time label aggregation.
/// First-seen order is stored explicitly so hash iteration order cannot affect
/// ties or the original deterministic behavior.
fn best_label<'a>(
    atoms: impl Iterator<Item = AtomRef<'a>>,
    scores: &mut HashMap<AltId, LabelScore>,
) -> LabelChoice {
    scores.clear();

    for atom in atoms {
        let Some(label) = atom.alt_id() else {
            continue;
        };

        if label.is_blank() {
            continue;
        }

        let occupancy = occupancy(atom);
        let order = scores.len();

        scores
            .entry(label)
            .and_modify(|score| score.occupancy += occupancy)
            .or_insert(LabelScore { occupancy, order });
    }

    LabelChoice {
        selected: best_scored_label(scores),
        has_altlocs: !scores.is_empty(),
    }
}

/// Selects the winning label from accumulated occupancy scores.
///
/// Selection reproduces ordered `>` semantics exactly: ties retain the first
/// label encountered, and a NaN accumulated by the first label continues to
/// dominate exactly as it did under the original sequential comparison.
fn best_scored_label(scores: &HashMap<AltId, LabelScore>) -> Option<AltId> {
    let mut earliest: Option<(AltId, LabelScore)> = None;
    let mut best_finite: Option<(AltId, LabelScore)> = None;

    for (&label, &score) in scores {
        if earliest.is_none_or(|(_, current)| score.order < current.order) {
            earliest = Some((label, score));
        }

        if score.occupancy.is_nan() {
            continue;
        }

        if best_finite.is_none_or(|(_, current)| {
            score.occupancy > current.occupancy
                || (score.occupancy.to_bits() == current.occupancy.to_bits()
                    && score.order < current.order)
        }) {
            best_finite = Some((label, score));
        }
    }

    match earliest {
        Some((label, score)) if score.occupancy.is_nan() => Some(label),
        Some(_) => best_finite.map(|(label, _)| label),
        None => None,
    }
}

/// Updates the retained alternate row for one atom name.
///
/// Equal occupancies preserve the first encountered atom because replacement
/// occurs only for a strictly greater occupancy.
fn record_atom_choice(
    choices: &mut HashMap<SymbolId, AtomChoice>,
    name: SymbolId,
    atom: AtomRef<'_>,
) {
    let candidate = AtomChoice {
        atom: atom.index().get(),
        occupancy: occupancy(atom),
    };

    choices
        .entry(name)
        .and_modify(|current| {
            if candidate.occupancy > current.occupancy {
                *current = candidate;
            }
        })
        .or_insert(candidate);
}

/// Returns the recorded atom occupancy, defaulting absent values to one.
fn occupancy(atom: AtomRef<'_>) -> f32 {
    match atom.occupancy() {
        Some(value) => value,
        None => 1.0,
    }
}

/// Returns whether an atom belongs to the selected conformation.
///
/// Blank-labelled atoms belong to every conformation and are always retained.
fn keeps_label(atom: AtomRef<'_>, selected: Option<AltId>) -> bool {
    atom.alt_id()
        .is_some_and(|label| label.is_blank() || Some(label) == selected)
}

/// Produces owned diagnostic context for a named alternate-location policy.
///
/// Allocation occurs only on the missing-label warning path, where the
/// diagnostic must retain its own copy of the requested identifier.
fn named_policy_label(policy: &AltlocPolicy) -> Box<str> {
    match policy {
        AltlocPolicy::Label(label) => label.clone(),
        _ => "".into(),
    }
}

#[cfg(test)]
#[path = "disorder_tests.rs"]
mod tests;
