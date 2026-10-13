//! Borrowed readers for component atom, bond and coordinate rows.
//!
//! A view holds its names as text borrowed from the category, so reading a row
//! allocates only when the file wrote an identifier as a number. Building a
//! [`Component`](crate::Component) copies the views once; checking that a row
//! is well formed copies nothing. Both go through the same readers, so a row
//! the check accepts is a row the lowering accepts. Reading a row costs time
//! proportional to its item count.

use super::items::{
    formal_charge, identifier_ref, invalid_value, one_finding, required_flag, required_number,
    text_ref, unknown_marker,
};
use crate::model::{BondStereo, StereoConfiguration};
use molframe_cif::Category;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{BondOrder, Element};
use num_traits::ToPrimitive;
use std::borrow::Cow;

type LoweringResult<T> = Result<T, Vec<Diagnostic>>;

/// One atom row read in place.
pub(super) struct AtomView<'a> {
    pub(super) name: Cow<'a, str>,
    pub(super) alternate_name: Option<Cow<'a, str>>,
    pub(super) element: Element,
    pub(super) charge: i8,
    pub(super) aromatic: bool,
    pub(super) leaving: bool,
    pub(super) stereo: Option<StereoConfiguration>,
}

/// One bond row read in place.
pub(super) struct BondView<'a> {
    pub(super) atom_a: Cow<'a, str>,
    pub(super) atom_b: Cow<'a, str>,
    pub(super) order: BondOrder,
    pub(super) aromatic: bool,
    pub(super) stereo: Option<BondStereo>,
}

/// Reads an atom row; `lenient` accepts the reduced rows of a deposited
/// structure file, which omit the formal charge and the leaving-atom flag.
pub(super) fn atom_view(
    category: &Category,
    row: usize,
    lenient: bool,
) -> Result<AtomView<'_>, Diagnostic> {
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
    Ok(AtomView {
        name: identifier_ref(category, "atom_id", row)?,
        alternate_name: category.identifier("alt_atom_id", row),
        element,
        charge,
        aromatic: required_flag(category, "pdbx_aromatic_flag", row)?,
        leaving,
        stereo: optional_stereo(category, row)?,
    })
}

/// Reads a bond row.
pub(super) fn bond_view(category: &Category, row: usize) -> Result<BondView<'_>, Diagnostic> {
    let aromatic = required_flag(category, "pdbx_aromatic_flag", row)?;
    Ok(BondView {
        atom_a: identifier_ref(category, "atom_id_1", row)?,
        atom_b: identifier_ref(category, "atom_id_2", row)?,
        order: bond_order(category, row, aromatic)?,
        aromatic,
        stereo: bond_stereo(category, row)?,
    })
}

/// The element of a component atom.
///
/// The archive writes `X` for an atom whose element is deliberately
/// unspecified, as in the ambiguous residues ASX and GLX. That is reported as
/// [`Element::UNKNOWN`] rather than guessed or rejected.
fn element(category: &Category, row: usize) -> Result<Element, Diagnostic> {
    let symbol = text_ref(category, "type_symbol", row)?;
    if symbol.trim().eq_ignore_ascii_case("X") {
        return Ok(Element::UNKNOWN);
    }
    Element::from_symbol(symbol).ok_or_else(|| invalid_value(category, "type_symbol", row))
}

/// One coordinate set (ideal or model) for every atom of a component, handed
/// to `sink` row by row.
///
/// The archive writes `?` for a coordinate it does not have, and the two sets
/// are curated independently: a component can carry ideal coordinates and no
/// model coordinates, or the reverse. A set with any unknown or inapplicable
/// coordinate is reported as absent as a whole, because a partial set cannot be
/// indexed by atom and substituting zeros would place atoms at the origin. The
/// caller chooses between the sets; nothing here falls back from one to the
/// other. Returns whether the whole set was present; positions already given
/// to `sink` are meaningless when it was not.
pub(super) fn coordinate_rows(
    category: &Category,
    rows: impl Iterator<Item = usize>,
    ideal: bool,
    mut sink: impl FnMut([f32; 3]),
) -> LoweringResult<bool> {
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
        return Ok(false);
    }
    for row in rows {
        let mut position = [0.0; 3];
        for (axis, item) in items.iter().enumerate() {
            if unknown_marker(category, item, row) {
                return Ok(false);
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
        sink(position);
    }
    Ok(true)
}

fn bond_order(category: &Category, row: usize, aromatic: bool) -> Result<BondOrder, Diagnostic> {
    if aromatic {
        return Ok(BondOrder::Aromatic);
    }
    let text = text_ref(category, "value_order", row)?;
    let is = |name: &str| text.eq_ignore_ascii_case(name);
    if is("SING") {
        Ok(BondOrder::Single)
    } else if is("DOUB") {
        Ok(BondOrder::Double)
    } else if is("TRIP") {
        Ok(BondOrder::Triple)
    } else if is("QUAD") {
        Ok(BondOrder::Quadruple)
    } else {
        Err(invalid_value(category, "value_order", row))
    }
}

/// Bond stereo: `N` is none, and `E`/`Z` are double-bond geometry.
///
/// Anything else is rejected rather than reinterpreted: the tetrahedral values
/// belong to atoms, and a bond carrying one is a malformed dictionary entry.
fn bond_stereo(category: &Category, row: usize) -> Result<Option<BondStereo>, Diagnostic> {
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
    category: &Category,
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

/// Whether a component's atoms and bonds form a valid graph: unique atom names,
/// bonds between known atoms, no self-bond and no repeated pair.
///
/// Sorts both slices in place, so the cost is `O((a + b) log (a + b))` with no
/// allocation.
pub(super) fn topology_is_valid(atoms: &mut [AtomView<'_>], bonds: &mut [BondView<'_>]) -> bool {
    atoms.sort_unstable_by(|left, right| left.name.cmp(&right.name));
    if atoms.windows(2).any(|pair| pair[0].name == pair[1].name) {
        return false;
    }
    let known = |name: &str| {
        atoms
            .binary_search_by(|atom| atom.name.as_ref().cmp(name))
            .is_ok()
    };
    if bonds
        .iter()
        .any(|bond| !known(&bond.atom_a) || !known(&bond.atom_b) || bond.atom_a == bond.atom_b)
    {
        return false;
    }
    bonds.sort_unstable_by(|left, right| ordered(left).cmp(&ordered(right)));
    !bonds
        .windows(2)
        .any(|pair| ordered(&pair[0]) == ordered(&pair[1]))
}

fn ordered<'a>(bond: &'a BondView<'_>) -> (&'a str, &'a str) {
    if bond.atom_a <= bond.atom_b {
        (&bond.atom_a, &bond.atom_b)
    } else {
        (&bond.atom_b, &bond.atom_a)
    }
}

/// Reports that a component's atom count disagrees with its coordinate count.
pub(super) fn coordinate_count_mismatch(
    category: &Category,
    expected: usize,
    observed: usize,
) -> Vec<Diagnostic> {
    vec![
        Diagnostic::new(Code::E3001)
            .in_category(category.name())
            .with_context("expected", expected.to_string())
            .with_context("observed", observed.to_string()),
    ]
}
