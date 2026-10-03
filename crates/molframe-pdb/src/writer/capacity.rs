//! What the fixed-column format cannot hold, checked before any byte is written.

use super::{PdbOptions, atom_name, chain_label, component_name, residue_sequence};
use crate::hybrid36;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::Select;
use molframe_core::structure::{AtomRef, ResidueRef, Structure};

/// Bounds whose rounding would overflow an eight-column, three-decimal field.
/// The negative side loses one digit to the sign, so the limits are asymmetric.
const MIN_COORDINATE: f64 = -999.999_5;
const MAX_COORDINATE: f64 = 9_999.999_5;

#[derive(Clone, Copy)]
pub(crate) enum RequiredAtomFields {
    Pdb,
    Pqr,
    Pdbqt,
}

/// Everything about the structure that the format cannot hold.
pub(crate) fn check_capacity(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
    required: RequiredAtomFields,
) -> Vec<Diagnostic> {
    let data = structure.data();
    let mut refusals = Vec::new();

    for chain in data.chains() {
        if !select.accept_chain(chain.index()) {
            continue;
        }
        let original = chain_label(&chain, options.identifier_namespace());
        let written = options.label_for(original);
        if written.chars().count() > 1 {
            refusals.push(
                Diagnostic::new(Code::E4102)
                    .with_context("chain", original.to_string())
                    .with_context("written as", written.to_string()),
            );
        }
        for residue in chain.residues() {
            if !select.accept_residue(residue.index()) {
                continue;
            }
            check_residue(&residue, options, select, required, &mut refusals);
        }
    }

    let atoms = i64::from(structure.atom_count());
    if !options.uses_hybrid36() && hybrid36::needs_hybrid36(atoms, 5) {
        refusals.push(
            Diagnostic::new(Code::E4101)
                .with_context("atoms", atoms.to_string())
                .with_context("field holds", "99999"),
        );
    }
    refusals
}

fn check_residue(
    residue: &ResidueRef<'_>,
    options: &PdbOptions,
    select: &impl Select,
    required: RequiredAtomFields,
    refusals: &mut Vec<Diagnostic>,
) {
    if let Some(seq) = residue_sequence(*residue, options.identifier_namespace())
        && !options.uses_hybrid36()
        && hybrid36::needs_hybrid36(i64::from(seq), 4)
    {
        refusals.push(
            Diagnostic::new(Code::E4103)
                .with_context("residue", seq.to_string())
                .with_context("field holds", "9999"),
        );
    }
    for atom in residue.atoms() {
        if !select.accept_atom(atom.index()) {
            continue;
        }
        check_required_atom_fields(atom, *residue, options, required, refusals);
        let Some(position) = atom.position() else {
            continue;
        };
        if position.iter().any(|value| {
            let value = f64::from(*value);
            !value.is_finite() || !(MIN_COORDINATE..MAX_COORDINATE).contains(&value)
        }) {
            refusals.push(
                Diagnostic::new(Code::E4104)
                    .with_context("atom", atom.index().to_string())
                    .with_context("position", format!("{position:?}")),
            );
        }
    }
}

fn check_required_atom_fields(
    atom: AtomRef<'_>,
    residue: ResidueRef<'_>,
    options: &PdbOptions,
    required: RequiredAtomFields,
    refusals: &mut Vec<Diagnostic>,
) {
    let common = [
        (atom.position().is_some(), "coordinates"),
        (
            atom_name(atom, options.identifier_namespace()).is_some(),
            "atom identifier",
        ),
        (
            component_name(atom, residue, options.identifier_namespace()).is_some(),
            "component identifier",
        ),
        (
            residue_sequence(residue, options.identifier_namespace()).is_some(),
            "residue sequence identifier",
        ),
    ];
    for (present, field) in common {
        if !present {
            refusals.push(missing_field(field, atom, options));
        }
    }
    if matches!(
        required,
        RequiredAtomFields::Pdb | RequiredAtomFields::Pdbqt
    ) {
        for (present, field) in [
            (atom.occupancy().is_some(), "occupancy"),
            (atom.b_factor().is_some(), "B factor"),
        ] {
            if !present {
                refusals.push(missing_field(field, atom, options));
            }
        }
    }
}

pub(super) fn required_text<'a>(
    value: Option<&'a str>,
    field: &'static str,
    atom: AtomRef<'_>,
    options: &PdbOptions,
) -> Result<&'a str, Diagnostic> {
    value.ok_or_else(|| missing_field(field, atom, options))
}

pub(super) fn required_value<T>(
    value: Option<T>,
    field: &'static str,
    atom: AtomRef<'_>,
    options: &PdbOptions,
) -> Result<T, Diagnostic> {
    value.ok_or_else(|| missing_field(field, atom, options))
}

fn missing_field(field: &'static str, atom: AtomRef<'_>, options: &PdbOptions) -> Diagnostic {
    Diagnostic::new(Code::E4105)
        .with_context("required field", field)
        .with_context("atom", atom.index().to_string())
        .with_context("namespace", options.identifier_namespace().as_str())
}
