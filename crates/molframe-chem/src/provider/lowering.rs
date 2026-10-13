//! Lowering of CCD data blocks into [`Component`] values.
//!
//! Each reader reports every finding it can before giving up, so a malformed
//! dictionary entry is diagnosed whole rather than one field at a time.

use super::items::{collect_rows, one_finding, required_identifier, required_text};
use super::validation::validate_component_topology;
use super::views::{atom_view, bond_view, coordinate_count_mismatch, coordinate_rows};
use crate::model::{Component, ComponentAtom, ComponentBond, ComponentKind};
use molframe_cif::{Category, DataBlock};
use molframe_core::Element;
use molframe_core::diagnostic::{Code, Diagnostic};
use std::sync::Arc;

type ComponentCoordinates = Arc<[[f32; 3]]>;
type LoweringResult<T> = Result<T, Vec<Diagnostic>>;

/// The rows of one category that belong to a single component.
///
/// A dictionary block holds one component, so its rows are every row; an
/// in-file dictionary interleaves components, so its rows are those whose
/// `comp_id` matches.
#[derive(Clone, Copy)]
pub(super) struct Rows<'a> {
    pub(super) category: &'a Category,
    pub(super) rows: &'a [usize],
}

pub(super) fn component(block: &DataBlock) -> LoweringResult<Component> {
    let Some(metadata) = block.category("chem_comp") else {
        return Err(vec![
            Diagnostic::new(Code::E2001)
                .with_context("block", block.name())
                .in_category("chem_comp"),
        ]);
    };
    let all = |name: &str| {
        block
            .category(name)
            .map(|category| (category, (0..category.row_count()).collect::<Vec<_>>()))
    };
    let atoms = all("chem_comp_atom");
    let bonds = all("chem_comp_bond");
    component_from(
        metadata,
        0,
        atoms
            .as_ref()
            .map(|(category, rows)| Rows { category, rows }),
        bonds
            .as_ref()
            .map(|(category, rows)| Rows { category, rows }),
        false,
    )
}

/// Lowers one component with its rows copied out of the category.
///
/// Time and allocation are proportional to the component's own rows.
pub(super) fn component_from(
    metadata: &Category,
    row: usize,
    atom_rows: Option<Rows<'_>>,
    bond_rows: Option<Rows<'_>>,
    lenient: bool,
) -> LoweringResult<Component> {
    let id = required_identifier(metadata, "id", row).map_err(one_finding)?;
    let name = required_identifier(metadata, "name", row).map_err(one_finding)?;
    let component_type = required_text(metadata, "type", row).map_err(one_finding)?;
    let atoms = atoms(atom_rows, lenient)?;
    let bonds = bonds(bond_rows)?;
    validate_component_topology(&atoms, &bonds)?;
    let component_kind = classify_component(&component_type, &atoms, &bonds);
    let ideal_coordinates = coordinates(atom_rows, true, atoms.len())?;
    let model_coordinates = coordinates(atom_rows, false, atoms.len())?;
    Ok(Component {
        id,
        name,
        kind: component_kind,
        parent: metadata
            .identifier("mon_nstd_parent_comp_id", row)
            .map(|value| value.into_owned().into()),
        one_letter_code: metadata
            .identifier("one_letter_code", row)
            .and_then(|value| one_letter_code(value.as_ref())),
        formula: metadata
            .identifier("formula", row)
            .map(|value| value.into_owned().into()),
        ideal_coordinates,
        model_coordinates,
        atoms: atoms.into(),
        bonds: bonds.into(),
    })
}

pub(super) fn classify_component(
    component_type: &str,
    atoms: &[ComponentAtom],
    bonds: &[ComponentBond],
) -> ComponentKind {
    let broad = kind(Some(component_type));
    if broad == ComponentKind::NonPolymer && is_water(atoms) {
        return ComponentKind::Solvent;
    }
    if broad == ComponentKind::NonPolymer
        && atoms.len() == 1
        && bonds.is_empty()
        && atoms[0].charge != 0
    {
        ComponentKind::Ion
    } else {
        broad
    }
}

/// A neutral molecule of one oxygen and two hydrogen (or deuterium) atoms.
///
/// The archive types water, heavy water and their relatives `NON-POLYMER`, so
/// the type alone would leave solvent indistinguishable from a ligand.
fn is_water(atoms: &[ComponentAtom]) -> bool {
    let oxygen = Element::from_symbol("O");
    let hydrogen = [Element::from_symbol("H"), Element::from_symbol("D")];
    atoms.len() == 3
        && atoms.iter().all(|atom| atom.charge == 0)
        && atoms
            .iter()
            .filter(|atom| Some(atom.element) == oxygen)
            .count()
            == 1
        && atoms
            .iter()
            .filter(|atom| hydrogen.contains(&Some(atom.element)))
            .count()
            == 2
}

fn one_letter_code(value: &str) -> Option<u8> {
    let bytes = value.as_bytes();
    (bytes.len() == 1 && bytes[0].is_ascii()).then_some(bytes[0])
}

pub(super) fn kind(value: Option<&str>) -> ComponentKind {
    const PEPTIDE: &str = "PEPTIDE";
    const DNA: &str = "DNA";
    const RNA: &str = "RNA";
    const SACCHARIDE: &str = "SACCHARIDE";
    const LIPID: &str = "LIPID";
    const WATER: &str = "WATER";
    const NON_POLYMER: &str = "NON-POLYMER";

    let Some(text) = value else {
        return ComponentKind::Unknown;
    };
    let normalized = text.trim().to_ascii_uppercase();
    let has_word = |word: &str| {
        normalized
            .split(|character: char| !character.is_ascii_alphanumeric())
            .any(|token| token == word)
    };
    if has_word(PEPTIDE) {
        ComponentKind::AminoAcid
    } else if has_word(DNA) || has_word(RNA) {
        ComponentKind::Nucleotide
    } else if has_word(SACCHARIDE) {
        ComponentKind::Saccharide
    } else if has_word(LIPID) {
        ComponentKind::Lipid
    } else if normalized == WATER {
        ComponentKind::Solvent
    } else if normalized == NON_POLYMER {
        ComponentKind::NonPolymer
    } else {
        ComponentKind::Unknown
    }
}

/// Atoms of a component; `lenient` accepts the reduced rows of a deposited
/// structure file, which omit the formal charge and the leaving-atom flag.
fn atoms(rows: Option<Rows<'_>>, lenient: bool) -> LoweringResult<Vec<ComponentAtom>> {
    // A component with no atom category is a legitimate entry with no atoms,
    // such as the unknown-ligand placeholder; an empty category is not.
    let Some(Rows { category, rows }) = rows else {
        return Ok(Vec::new());
    };
    if rows.is_empty() {
        return Err(vec![
            Diagnostic::new(Code::E2002)
                .with_context("item", "atom_id")
                .with_context("row", "0")
                .in_category(category.name()),
        ]);
    }
    collect_rows(rows, |row| {
        let view = atom_view(category, row, lenient)?;
        Ok(ComponentAtom {
            name: view.name.into_owned().into(),
            alternate_name: view.alternate_name.map(|name| name.into_owned().into()),
            element: view.element,
            charge: view.charge,
            aromatic: view.aromatic,
            leaving: view.leaving,
            stereo: view.stereo,
        })
    })
}

fn bonds(rows: Option<Rows<'_>>) -> LoweringResult<Vec<ComponentBond>> {
    let Some(Rows { category, rows }) = rows else {
        return Ok(Vec::new());
    };
    collect_rows(rows, |row| {
        let view = bond_view(category, row)?;
        Ok(ComponentBond {
            atom_a: view.atom_a.into_owned().into(),
            atom_b: view.atom_b.into_owned().into(),
            order: view.order,
            aromatic: view.aromatic,
            stereo: view.stereo,
        })
    })
}

fn coordinates(
    atoms: Option<Rows<'_>>,
    ideal: bool,
    expected: usize,
) -> LoweringResult<Option<ComponentCoordinates>> {
    let Some(Rows { category, rows }) = atoms else {
        return Ok(None);
    };
    let mut positions = Vec::with_capacity(expected);
    let present = coordinate_rows(category, rows.iter().copied(), ideal, |position| {
        positions.push(position);
    })?;
    if !present {
        return Ok(None);
    }
    if positions.len() == expected {
        Ok(Some(positions.into()))
    } else {
        Err(coordinate_count_mismatch(
            category,
            expected,
            positions.len(),
        ))
    }
}
