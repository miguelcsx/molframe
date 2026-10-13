//! Fixed-column fields shared by the atom, `ANISOU`, `LINK` and `TER` records.

use super::pdb::{PdbOptions, component_name, insertion_code, residue_field, residue_sequence};
use molframe_core::structure::{AtomRef, ResidueRef};

/// The four-column name field.
///
/// A name shorter than four characters starts in the second column, which is the
/// convention that lets a reader tell a one-letter element from a two-letter one.
pub(crate) fn atom_name_field(name: &str, element: &str) -> String {
    let two_letter = element.len() == 2
        && name
            .get(..2)
            .is_some_and(|head| head.eq_ignore_ascii_case(element));
    if name.len() >= 4 {
        name.to_owned()
    } else if two_letter {
        format!("{name:<4}")
    } else {
        format!(" {name:<3}")
    }
}

/// The element symbol in the upper case columns 77-78 use.
pub(crate) fn element_column(atom: &AtomRef<'_>) -> String {
    atom.element()
        .map_or_else(String::new, |element| element.symbol().to_ascii_uppercase())
}

/// Columns 18-27 of `TER`: the last polymer residue the chain ended on.
pub(super) fn terminus_fields(
    residue: &ResidueRef<'_>,
    chain: &str,
    options: &PdbOptions,
) -> Option<String> {
    let namespace = options.identifier_namespace();
    let atom = residue.atoms().next()?;
    let component = component_name(atom, *residue, namespace)?;
    let seq = i64::from(residue_sequence(*residue, namespace)?);
    let fields = format!(
        "      {component:>3} {chain:>1}{}{}",
        residue_field(seq, options),
        insertion_code(*residue)
    );
    Some(fields.trim_end().to_owned())
}

#[cfg(test)]
#[path = "columns_tests.rs"]
mod tests;
