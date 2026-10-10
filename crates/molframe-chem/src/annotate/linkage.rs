//! Inter-residue polymer bond inference.
//!
//! Every eligible attachment rule is evaluated, atom by atom, so a rule list
//! never stops at the first rule whose component kinds match. Candidate atom
//! pairs must exist, be alternate-location compatible, and sit inside a
//! bond-type-specific distance window. Competing candidates are then resolved
//! per conformer state: the shortest candidate wins each atom, and a genuine
//! tie is reported instead of being guessed.

use super::{PolymerLinkPolicy, PolymerLinkRule, alt_compatible};
use crate::Component;
use molframe_core::BondOrder;
use molframe_core::Element;
use molframe_core::bond::{BondProvenance, BondRecord, BondTableBuilder};
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::structure::{AtomRef, ResidueRef};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Distances closer than this (angstrom) cannot be told apart.
const TIE_TOLERANCE: f32 = 0.05;
/// No covalent bond between heavy atoms is shorter than this (angstrom).
const GENERIC_MIN: f32 = 0.8;

/// The inclusive covalent-bond distance window for a bonded element pair.
///
/// Known pairs use a tight window so a clash or a copy-paste coordinate error
/// is not mistaken for a bond. Unknown pairs only get a generic lower bound and
/// are limited from above by the policy distance.
fn bond_window(a: Option<Element>, b: Option<Element>) -> (f32, f32) {
    let pair =
        |x: Element, y: Element| (a == Some(x) && b == Some(y)) || (a == Some(y) && b == Some(x));
    if pair(Element::CARBON, Element::NITROGEN) {
        (1.2, 1.6)
    } else if pair(Element::CARBON, Element::OXYGEN) {
        (1.3, 1.6)
    } else if pair(Element::SULFUR, Element::SULFUR) {
        (1.8, 2.3)
    } else if pair(Element::PHOSPHORUS, Element::OXYGEN) {
        (1.4, 1.8)
    } else {
        (GENERIC_MIN, f32::INFINITY)
    }
}

fn distance(a: AtomRef<'_>, b: AtomRef<'_>) -> Option<f32> {
    let (a, b) = (a.position()?, b.position()?);
    let squared: f32 = a.iter().zip(b).map(|(l, r)| (l - r).powi(2)).sum();
    Some(squared.sqrt())
}

/// The shared conformer state of a pair: the non-blank label either atom has.
fn state_of(a: AtomRef<'_>, b: AtomRef<'_>) -> Option<molframe_core::AltId> {
    [a.alt_id(), b.alt_id()]
        .into_iter()
        .flatten()
        .find(|alt| !alt.is_blank())
}

struct Candidate<'a> {
    left: AtomRef<'a>,
    right: AtomRef<'a>,
    distance: f32,
    state: Option<molframe_core::AltId>,
}

fn shares_atom(a: &Candidate<'_>, b: &Candidate<'_>) -> bool {
    a.left.index() == b.left.index() || a.right.index() == b.right.index()
}

/// Whether the candidate lies within the window, recording why when it does not.
enum Verdict {
    Valid,
    TooShort,
    TooLong,
}

fn judge(distance: f32, atom_a: AtomRef<'_>, atom_b: AtomRef<'_>, max: f32) -> Verdict {
    let (low, high) = bond_window(atom_a.element(), atom_b.element());
    if distance < low {
        Verdict::TooShort
    } else if distance > high.min(max) {
        Verdict::TooLong
    } else {
        Verdict::Valid
    }
}

fn named_atoms<'a>(residue: ResidueRef<'a>, name: &str) -> Vec<AtomRef<'a>> {
    residue
        .atoms()
        .filter(|atom| atom.name() == Some(name))
        .collect()
}

pub(super) fn link_polymer(
    residues: &[ResidueRef<'_>],
    components: &BTreeMap<Box<str>, Arc<Component>>,
    bonds: &mut BondTableBuilder,
    findings: &mut Diagnostics,
    policy: &PolymerLinkPolicy,
) {
    let PolymerLinkPolicy::Explicit { angstrom, rules } = policy else {
        return;
    };
    for pair in residues.windows(2) {
        let [left, right] = pair else {
            continue;
        };
        let Some(left_component) = left.name().and_then(|name| components.get(name)) else {
            continue;
        };
        let Some(right_component) = right.name().and_then(|name| components.get(name)) else {
            continue;
        };
        let eligible: Vec<&PolymerLinkRule> = rules
            .iter()
            .filter(|rule| {
                rule.left_kind == left_component.kind && rule.right_kind == right_component.kind
            })
            .collect();
        if eligible.is_empty() {
            continue;
        }
        link_pair(*left, *right, &eligible, *angstrom, bonds, findings);
    }
}

fn link_pair(
    left: ResidueRef<'_>,
    right: ResidueRef<'_>,
    rules: &[&PolymerLinkRule],
    max: f32,
    bonds: &mut BondTableBuilder,
    findings: &mut Diagnostics,
) {
    let mut valid: Vec<Candidate<'_>> = Vec::new();
    let (mut examined, mut too_short) = (false, false);
    for rule in rules {
        for atom_a in named_atoms(left, &rule.left_atom) {
            for atom_b in named_atoms(right, &rule.right_atom) {
                if !alt_compatible(atom_a, atom_b) {
                    continue;
                }
                let Some(d) = distance(atom_a, atom_b) else {
                    continue;
                };
                examined = true;
                match judge(d, atom_a, atom_b, max) {
                    Verdict::Valid => {
                        let duplicate = valid.iter().any(|c| {
                            c.left.index() == atom_a.index() && c.right.index() == atom_b.index()
                        });
                        if !duplicate {
                            valid.push(Candidate {
                                left: atom_a,
                                right: atom_b,
                                distance: d,
                                state: state_of(atom_a, atom_b),
                            });
                        }
                    }
                    Verdict::TooShort => too_short = true,
                    Verdict::TooLong => {}
                }
            }
        }
    }
    if valid.is_empty() {
        if examined {
            findings.push(
                Diagnostic::new(Code::W3302)
                    .with_context("left residue", left.index().to_string())
                    .with_context("right residue", right.index().to_string())
                    .with_context("reason", if too_short { "too short" } else { "too far" }),
            );
        }
        return;
    }
    let mut ambiguous = false;
    // Candidates only compete inside one conformer state: A-A and B-B are two
    // separate bonds, never one bond to be chosen between.
    let keep: Vec<bool> = valid
        .iter()
        .map(|c| {
            let rivals = valid.iter().filter(|o| {
                o.state == c.state && shares_atom(c, o) && o.distance + TIE_TOLERANCE < c.distance
            });
            rivals.count() == 0
        })
        .collect();
    let survivors: Vec<&Candidate<'_>> = valid
        .iter()
        .zip(&keep)
        .filter_map(|(c, keep)| keep.then_some(c))
        .collect();
    let mut tied = vec![false; survivors.len()];
    for (i, a) in survivors.iter().enumerate() {
        for (j, b) in survivors.iter().enumerate() {
            if i < j && a.state == b.state && shares_atom(a, b) {
                ambiguous = true;
                mark(&mut tied, i);
                mark(&mut tied, j);
            }
        }
    }
    for (candidate, tied) in survivors.iter().zip(&tied) {
        if *tied {
            continue;
        }
        bonds.push(BondRecord {
            atom_a: candidate.left.index(),
            atom_b: candidate.right.index(),
            order: BondOrder::Polymeric,
            provenance: BondProvenance::ChemicalComponentDictionary,
        });
    }
    if ambiguous {
        findings.push(
            Diagnostic::new(Code::W3302)
                .with_context("left residue", left.index().to_string())
                .with_context("right residue", right.index().to_string())
                .with_context("reason", "ambiguous attachment candidates"),
        );
    }
}

fn mark(flags: &mut [bool], at: usize) {
    if let Some(flag) = flags.get_mut(at) {
        *flag = true;
    }
}
