//! Lowering of CCD data blocks into [`Component`] values.
//!
//! Each reader reports every finding it can before giving up, so a malformed
//! dictionary entry is diagnosed whole rather than one field at a time.

use super::items::{
    collect_rows, formal_charge, invalid_value, one_finding, required_flag, required_identifier,
    required_number, required_text, unknown_marker,
};
use super::validation::validate_component_topology;
use crate::model::{
    BondStereo, Component, ComponentAtom, ComponentBond, ComponentKind, StereoConfiguration,
};
use molframe_cif::{Category, DataBlock};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{BondOrder, Element};
use num_traits::ToPrimitive;
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
        let element = element(category, row)?;
        let charge = if lenient && category.column("charge").is_none() {
            0
        } else {
            formal_charge(category, row)?
        };
        let leaving = if lenient && category.column("pdbx_leaving_atom_flag").is_none() {
            false
        } else {
            required_flag(category, "pdbx_leaving_atom_flag", row)?
        };
        Ok(ComponentAtom {
            name: required_identifier(category, "atom_id", row)?,
            alternate_name: category
                .identifier("alt_atom_id", row)
                .map(|value| value.into_owned().into()),
            element,
            charge,
            aromatic: required_flag(category, "pdbx_aromatic_flag", row)?,
            leaving,
            stereo: optional_stereo(category, row)?,
        })
    })
}

/// The element of a component atom.
///
/// The archive writes `X` for an atom whose element is deliberately
/// unspecified, as in the ambiguous residues ASX and GLX. That is reported as
/// [`Element::UNKNOWN`] rather than guessed or rejected.
fn element(category: &molframe_cif::Category, row: usize) -> Result<Element, Diagnostic> {
    let symbol = required_text(category, "type_symbol", row)?;
    if symbol.trim().eq_ignore_ascii_case("X") {
        return Ok(Element::UNKNOWN);
    }
    Element::from_symbol(&symbol).ok_or_else(|| invalid_value(category, "type_symbol", row))
}

fn bonds(rows: Option<Rows<'_>>) -> LoweringResult<Vec<ComponentBond>> {
    let Some(Rows { category, rows }) = rows else {
        return Ok(Vec::new());
    };
    collect_rows(rows, |row| {
        let aromatic = required_flag(category, "pdbx_aromatic_flag", row)?;
        Ok(ComponentBond {
            atom_a: required_identifier(category, "atom_id_1", row)?,
            atom_b: required_identifier(category, "atom_id_2", row)?,
            order: bond_order(category, row, aromatic)?,
            aromatic,
            stereo: bond_stereo(category, row)?,
        })
    })
}

/// One coordinate set (ideal or model) for every atom of a component.
///
/// The archive writes `?` for a coordinate it does not have, and the two sets
/// are curated independently: a component can carry ideal coordinates and no
/// model coordinates, or the reverse. A set with any unknown or inapplicable
/// coordinate is reported as absent as a whole, because a partial set cannot be
/// indexed by atom and substituting zeros would place atoms at the origin. The
/// caller chooses between the sets; nothing here falls back from one to the
/// other.
fn coordinates(
    atoms: Option<Rows<'_>>,
    ideal: bool,
    expected: usize,
) -> LoweringResult<Option<ComponentCoordinates>> {
    let Some(Rows { category, rows }) = atoms else {
        return Ok(None);
    };
    let items = if ideal {
        [
            "pdbx_model_Cartn_x_ideal",
            "pdbx_model_Cartn_y_ideal",
            "pdbx_model_Cartn_z_ideal",
        ]
    } else {
        ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"]
    };
    if items.iter().all(|item| category.column(item).is_none()) {
        return Ok(None);
    }
    let mut positions = Vec::with_capacity(expected);
    for &row in rows {
        let mut position = [0.0; 3];
        for (axis, item) in items.iter().enumerate() {
            if unknown_marker(category, item, row) {
                return Ok(None);
            }
            let value = required_number(category, item, row).map_err(one_finding)?;
            if value < f64::from(f32::MIN) || value > f64::from(f32::MAX) {
                return Err(vec![invalid_value(category, item, row)]);
            }
            let Some(value) = value.to_f32() else {
                return Err(vec![invalid_value(category, item, row)]);
            };
            position[axis] = value;
        }
        positions.push(position);
    }
    if positions.len() == expected {
        Ok(Some(positions.into()))
    } else {
        Err(vec![
            Diagnostic::new(Code::E3001)
                .in_category(category.name())
                .with_context("expected", expected.to_string())
                .with_context("observed", positions.len().to_string()),
        ])
    }
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

fn bond_order(
    category: &molframe_cif::Category,
    row: usize,
    aromatic: bool,
) -> Result<BondOrder, Diagnostic> {
    if aromatic {
        return Ok(BondOrder::Aromatic);
    }
    match required_text(category, "value_order", row)?
        .to_ascii_uppercase()
        .as_str()
    {
        "SING" => Ok(BondOrder::Single),
        "DOUB" => Ok(BondOrder::Double),
        "TRIP" => Ok(BondOrder::Triple),
        "QUAD" => Ok(BondOrder::Quadruple),
        _ => Err(invalid_value(category, "value_order", row)),
    }
}

/// Bond stereo: `N` is none, and `E`/`Z` are double-bond geometry.
///
/// Anything else is rejected rather than reinterpreted: the tetrahedral values
/// belong to atoms, and a bond carrying one is a malformed dictionary entry.
fn bond_stereo(
    category: &molframe_cif::Category,
    row: usize,
) -> Result<Option<BondStereo>, Diagnostic> {
    let item = "pdbx_stereo_config";
    match category.text(item, row) {
        None => Ok(None),
        Some(value) if value.eq_ignore_ascii_case("N") => Ok(None),
        Some(value) if value.eq_ignore_ascii_case("E") => Ok(Some(BondStereo::E)),
        Some(value) if value.eq_ignore_ascii_case("Z") => Ok(Some(BondStereo::Z)),
        Some(_) => Err(invalid_value(category, item, row)),
    }
}

fn optional_stereo(
    category: &molframe_cif::Category,
    row: usize,
) -> Result<Option<StereoConfiguration>, Diagnostic> {
    let item = "pdbx_stereo_config";
    match category.text(item, row) {
        None => Ok(None),
        Some(value) if value.eq_ignore_ascii_case("N") => Ok(None),
        Some(value) if value.eq_ignore_ascii_case("R") => Ok(Some(StereoConfiguration::R)),
        Some(value) if value.eq_ignore_ascii_case("S") => Ok(Some(StereoConfiguration::S)),
        Some(value) if value.eq_ignore_ascii_case("MIXED") => Ok(Some(StereoConfiguration::Mixed)),
        Some(_) => Err(invalid_value(category, item, row)),
    }
}
