//! Selection-aware ranges; O(residues log residues) identity preflight, O(residues) output.

use super::pdb::{PdbOptions, chain_label, insertion_code, residue_field, residue_sequence};
use molframe_core::SecondaryStructure as Ss;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::Select;
use molframe_core::structure::{ChainRef, ResidueRef, Structure};
use std::fmt;

pub(super) fn preflight(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
) -> Vec<Diagnostic> {
    let mut count = 0;
    let mut refusals = Vec::new();
    ranges(structure, select, |_, first, last, _| {
        count += 1;
        if last.index().get() - first.index().get() + 1 > 99_999 {
            refusals.push(
                Diagnostic::new(Code::E2001).with_context("field", "HELIX length exceeds 99999"),
            );
        }
    });
    if count > 999 {
        refusals.push(
            Diagnostic::new(Code::E2001)
                .with_context("field", "secondary record serial exceeds 999"),
        );
    }
    if let Some(models) = structure.ragged_models()
        && models.iter().any(|model| {
            model
                .secondary_structure()
                .iter()
                .any(|state| state.is_helix() || state.is_strand())
        })
    {
        refusals.push(
            Diagnostic::new(Code::E2001)
                .with_context("field", "secondary records require a single model snapshot"),
        );
    }
    if count == 0 {
        return refusals;
    }
    let namespace = options.identifier_namespace();
    let maximum = if options.uses_hybrid36() {
        crate::hybrid36::decode("zzzz", 4)
    } else {
        Some(9999)
    };
    let mut identities = Vec::with_capacity(structure.data().topology.residues.len());
    for chain in structure
        .data()
        .chains()
        .filter(|chain| select.accept_chain(chain.index()))
    {
        let original = match super::pdb::chain_name(chain, namespace) {
            Some(label) => label,
            None => "",
        };
        let label = options.label_for(original);
        for residue in chain.residues().filter(|residue| {
            select.accept_residue(residue.index())
                && residue.atoms().any(|atom| select.accept_atom(atom.index()))
        }) {
            if let Some(sequence) = residue_sequence(residue, namespace) {
                identities.push((label, sequence, insertion_code(residue)));
            }
            let kind = match structure
                .secondary_structure()
                .get(residue.index().as_usize())
            {
                Some(kind) => *kind,
                None => Ss::Unknown,
            };
            if !kind.is_helix() && !kind.is_strand() {
                continue;
            }
            let name = match namespace {
                super::pdb::PdbIdentifierNamespace::Label => residue.name(),
                super::pdb::PdbIdentifierNamespace::Auth => residue.auth_name(),
            };
            let valid_name =
                name.is_some_and(|name| !name.is_empty() && name.len() <= 3 && name.is_ascii());
            let valid_sequence = residue_sequence(residue, namespace).is_some_and(|seq| {
                maximum.is_some_and(|maximum| (-999..=maximum).contains(&i64::from(seq)))
            });
            let insertion = insertion_code(residue);
            if !valid_name
                || !valid_sequence
                || super::pdb::chain_name(chain, namespace).is_none()
                || label.len() > 1
                || !label.is_ascii()
                || insertion.len() != 1
                || !insertion.is_ascii()
            {
                refusals.push(
                    Diagnostic::new(Code::E2001)
                        .with_context("field", "secondary endpoint identity does not fit PDB"),
                );
            }
        }
    }
    identities.sort_unstable();
    if identities.windows(2).any(|pair| pair[0] == pair[1]) {
        refusals.push(
            Diagnostic::new(Code::E2001)
                .with_context("field", "ambiguous secondary residue identity"),
        );
    }
    refusals
}

pub(super) fn write(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
) {
    let mut id = 0;
    ranges(structure, select, |chain, first, last, kind| {
        let namespace = options.identifier_namespace();
        let (Some(begin), Some(end)) = (
            residue_sequence(first, namespace),
            residue_sequence(last, namespace),
        ) else {
            return;
        };
        let (first_name, last_name) = match namespace {
            super::pdb::PdbIdentifierNamespace::Label => (first.name(), last.name()),
            super::pdb::PdbIdentifierNamespace::Auth => (first.auth_name(), last.auth_name()),
        };
        let (Some(first_name), Some(last_name)) = (first_name, last_name) else {
            return;
        };
        let label = options.label_for(chain_label(&chain, namespace));
        let begin = residue_field(i64::from(begin), options);
        let end = residue_field(i64::from(end), options);
        let first_code = insertion_code(first);
        let last_code = insertion_code(last);
        id += 1;
        if kind == Ss::Strand {
            let _ = writeln!(
                out,
                "SHEET  {id:>3} {id:>3} {:>2} {first_name:>3}{label}{begin}{first_code} {last_name:>3} {label}{end}{last_code}{:>2}",
                1, 0
            );
        } else {
            let class = match kind {
                Ss::AlphaHelix => 1,
                Ss::PiHelix => 3,
                Ss::ThreeTenHelix => 5,
                Ss::PolyProline => 10,
                _ => 2,
            };
            let length = last.index().get() - first.index().get() + 1;
            let _ = writeln!(
                out,
                "HELIX  {id:>3} {id:>3} {first_name:>3} {label} {begin}{first_code} {last_name:>3} {label} {end}{last_code}{class:>2}{:31}{length:>5}",
                ""
            );
        }
    });
}

fn ranges<'a>(
    structure: &'a Structure,
    select: &impl Select,
    mut visit: impl FnMut(ChainRef<'a>, ResidueRef<'a>, ResidueRef<'a>, Ss),
) {
    for chain in structure.data().chains() {
        if !select.accept_chain(chain.index()) {
            continue;
        }
        let mut pending: Option<(ResidueRef<'_>, ResidueRef<'_>, Ss)> = None;
        for residue in chain.residues() {
            let kind = match structure
                .secondary_structure()
                .get(residue.index().as_usize())
            {
                Some(state) => *state,
                None => Ss::Unknown,
            };
            let accepted = select.accept_residue(residue.index())
                && residue.atoms().any(|atom| select.accept_atom(atom.index()))
                && (kind.is_helix() || kind.is_strand());
            if let Some((first, last, old)) = pending {
                if accepted && old == kind {
                    pending = Some((first, residue, kind));
                    continue;
                }
                visit(chain, first, last, old);
            }
            pending = accepted.then_some((residue, residue, kind));
        }
        if let Some((first, last, kind)) = pending {
            visit(chain, first, last, kind);
        }
    }
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
