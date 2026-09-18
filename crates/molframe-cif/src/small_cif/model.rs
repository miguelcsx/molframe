//! Core-CIF small-molecule crystal structures.

use crate::{CifValue, DataBlock, Document};
use molframe_core::structure::UnitCell;
use std::collections::BTreeMap;

/// One crystallographic atom site from core CIF.
#[derive(Clone, Debug, PartialEq)]
pub struct SmallCifAtom {
    /// Unique atom-site label.
    pub label: Box<str>,
    /// Chemical type symbol as deposited.
    pub type_symbol: Box<str>,
    /// Fractional coordinates, when recorded.
    pub fractional: Option<[f64; 3]>,
    /// Cartesian coordinates in ångströms, when recorded.
    pub cartesian: Option<[f64; 3]>,
    /// Site occupancy.
    pub occupancy: Option<f64>,
    /// Isotropic displacement parameter U.
    pub u_iso: Option<f64>,
    /// Formal oxidation/charge text when present.
    pub charge: Option<Box<str>>,
}

/// One geometric bond record with resolved atom indices.
#[derive(Clone, Debug, PartialEq)]
pub struct SmallCifBond {
    /// First atom index.
    pub first: usize,
    /// Second atom index.
    pub second: usize,
    /// Reported bond distance in ångströms.
    pub distance: Option<f64>,
    /// Site-symmetry code applied to the second site.
    pub site_symmetry: Option<Box<str>>,
    /// Published bond/type marker.
    pub bond_type: Option<Box<str>>,
}

/// One core-CIF crystal model.
#[derive(Clone, Debug, PartialEq)]
pub struct SmallCifStructure {
    /// Data-block identifier.
    pub name: Box<str>,
    /// Unit cell.
    pub cell: Option<UnitCell>,
    /// Hermann-Mauguin space-group name.
    pub space_group_name: Option<Box<str>>,
    /// International Tables space-group number.
    pub space_group_number: Option<i32>,
    /// Explicit symmetry operations in source order.
    pub symmetry_operations: Vec<Box<str>>,
    /// Atom sites in source order.
    pub atoms: Vec<SmallCifAtom>,
    /// Resolved geometric bonds.
    pub bonds: Vec<SmallCifBond>,
}

/// Invalid core-CIF small-molecule model.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum SmallCifError {
    /// No data block or atom-site list is present.
    #[error("small-molecule CIF has no atom-site table")]
    MissingAtoms,
    /// An atom label is absent or duplicated.
    #[error("small-molecule CIF atom labels are absent or duplicated")]
    AtomLabels,
    /// Coordinate columns are incomplete or malformed.
    #[error("small-molecule CIF coordinates are incomplete or malformed")]
    Coordinates,
    /// Cell metadata is incomplete or invalid.
    #[error("small-molecule CIF unit cell is incomplete or invalid")]
    Cell,
    /// A bond references an absent atom label.
    #[error("small-molecule CIF bond references an absent atom")]
    BondAtom,
    /// A required deposited chemical type symbol is absent.
    #[error("small-molecule CIF atom type symbol is absent")]
    AtomTypeSymbol,
}

/// Dictionary dialect used to resolve core-CIF item names.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SmallCifDialect {
    /// Current DDL2 category/item names.
    #[default]
    Ddl2,
    /// Legacy DDL1 underscore tags, enabled only by explicit request.
    Ddl1,
}

/// Explicit decisions for lowering a small-molecule CIF.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SmallCifOptions {
    dialect: SmallCifDialect,
}

impl SmallCifOptions {
    /// Uses current DDL2 category/item names.
    #[must_use]
    pub const fn ddl2() -> Self {
        Self {
            dialect: SmallCifDialect::Ddl2,
        }
    }

    /// Uses legacy DDL1 underscore tags explicitly.
    #[must_use]
    pub const fn ddl1() -> Self {
        Self {
            dialect: SmallCifDialect::Ddl1,
        }
    }
}

/// Lowers the first core-CIF data block to a small-molecule crystal model.
/// Current DDL2 category/item names are used. Legacy DDL1 tags require
/// [`lower_small_cif_with_options`] and [`SmallCifOptions::ddl1`].
///
/// # Errors
///
/// Returns an error for missing/duplicate atom labels, incomplete coordinates
/// or cell data, or unresolved geometric bonds.
pub fn lower_small_cif(document: &Document) -> Result<SmallCifStructure, SmallCifError> {
    lower_small_cif_with_options(document, SmallCifOptions::ddl2())
}

/// Lowers the first core-CIF block under an explicit dictionary dialect.
///
/// # Errors
///
/// Returns an error for missing chemistry, identifiers, coordinates, cell data,
/// or bond endpoints under the selected dictionary dialect.
pub fn lower_small_cif_with_options(
    document: &Document,
    options: SmallCifOptions,
) -> Result<SmallCifStructure, SmallCifError> {
    let block = document.first_block().ok_or(SmallCifError::MissingAtoms)?;
    let source = Source {
        block,
        dialect: options.dialect,
    };
    let labels = source
        .column("atom_site", "label", "atom_site_label")
        .ok_or(SmallCifError::MissingAtoms)?;
    if labels.is_empty() {
        return Err(SmallCifError::MissingAtoms);
    }
    let mut atoms = Vec::with_capacity(labels.len());
    let mut label_indices = BTreeMap::new();
    for row in 0..labels.len() {
        let label = labels
            .get(row)
            .and_then(CifValue::as_identifier)
            .ok_or(SmallCifError::AtomLabels)?
            .into_owned();
        if label_indices.insert(label.clone(), row).is_some() {
            return Err(SmallCifError::AtomLabels);
        }
        let type_symbol = source
            .identifier("atom_site", "type_symbol", "atom_site_type_symbol", row)
            .ok_or(SmallCifError::AtomTypeSymbol)?
            .into();
        let fractional = coordinate(&source, "fract", row)?;
        let cartesian = coordinate(&source, "Cartn", row)?;
        if fractional.is_none() && cartesian.is_none() {
            return Err(SmallCifError::Coordinates);
        }
        atoms.push(SmallCifAtom {
            label: label.into(),
            type_symbol,
            fractional,
            cartesian,
            occupancy: source.number("atom_site", "occupancy", "atom_site_occupancy", row),
            u_iso: source.number(
                "atom_site",
                "U_iso_or_equiv",
                "atom_site_U_iso_or_equiv",
                row,
            ),
            charge: source
                .identifier("atom_site", "charge", "atom_site_charge", row)
                .map(Into::into),
        });
    }
    let cell = lower_cell(&source)?;
    let bonds = lower_bonds(&source, &label_indices)?;
    Ok(SmallCifStructure {
        name: block.name().into(),
        cell,
        space_group_name: match options.dialect {
            SmallCifDialect::Ddl2 => {
                source.identifier("space_group", "name_H-M_alt", "space_group_name_H-M_alt", 0)
            }
            SmallCifDialect::Ddl1 => source.identifier(
                "symmetry",
                "space_group_name_H-M",
                "symmetry_space_group_name_H-M",
                0,
            ),
        }
        .map(Into::into),
        space_group_number: match options.dialect {
            SmallCifDialect::Ddl2 => {
                source.integer("space_group", "IT_number", "space_group_IT_number", 0)
            }
            SmallCifDialect::Ddl1 => source.integer(
                "symmetry",
                "Int_Tables_number",
                "symmetry_Int_Tables_number",
                0,
            ),
        },
        symmetry_operations: symmetry_operations(&source),
        atoms,
        bonds,
    })
}

fn coordinate(
    source: &Source<'_>,
    namespace: &str,
    row: usize,
) -> Result<Option<[f64; 3]>, SmallCifError> {
    let values = ["x", "y", "z"].map(|axis| {
        let item = format!("{namespace}_{axis}");
        let tag = format!("atom_site_{namespace}_{axis}");
        source.number("atom_site", &item, &tag, row)
    });
    if values.iter().all(Option::is_none) {
        return Ok(None);
    }
    let [Some(x), Some(y), Some(z)] = values else {
        return Err(SmallCifError::Coordinates);
    };
    if [x, y, z].iter().any(|value| !value.is_finite()) {
        return Err(SmallCifError::Coordinates);
    }
    Ok(Some([x, y, z]))
}

fn lower_cell(source: &Source<'_>) -> Result<Option<UnitCell>, SmallCifError> {
    let values = [
        ("length_a", "cell_length_a"),
        ("length_b", "cell_length_b"),
        ("length_c", "cell_length_c"),
        ("angle_alpha", "cell_angle_alpha"),
        ("angle_beta", "cell_angle_beta"),
        ("angle_gamma", "cell_angle_gamma"),
    ]
    .map(|(item, tag)| source.number("cell", item, tag, 0));
    if values.iter().all(Option::is_none) {
        return Ok(None);
    }
    let [
        Some(a),
        Some(b),
        Some(c),
        Some(alpha),
        Some(beta),
        Some(gamma),
    ] = values
    else {
        return Err(SmallCifError::Cell);
    };
    if [a, b, c]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        || [alpha, beta, gamma]
            .iter()
            .any(|value| !value.is_finite() || !(0.0..180.0).contains(value))
    {
        return Err(SmallCifError::Cell);
    }
    Ok(Some(UnitCell {
        lengths: [a, b, c],
        angles: [alpha, beta, gamma],
    }))
}

fn lower_bonds(
    source: &Source<'_>,
    indices: &BTreeMap<String, usize>,
) -> Result<Vec<SmallCifBond>, SmallCifError> {
    let Some(first_column) = source.column(
        "geom_bond",
        "atom_site_label_1",
        "geom_bond_atom_site_label_1",
    ) else {
        return Ok(Vec::new());
    };
    let mut bonds = Vec::with_capacity(first_column.len());
    for row in 0..first_column.len() {
        let first_label = first_column
            .get(row)
            .and_then(CifValue::as_identifier)
            .ok_or(SmallCifError::BondAtom)?;
        let second_label = source
            .identifier(
                "geom_bond",
                "atom_site_label_2",
                "geom_bond_atom_site_label_2",
                row,
            )
            .ok_or(SmallCifError::BondAtom)?;
        bonds.push(SmallCifBond {
            first: *indices
                .get(first_label.as_ref())
                .ok_or(SmallCifError::BondAtom)?,
            second: *indices
                .get(second_label.as_ref())
                .ok_or(SmallCifError::BondAtom)?,
            distance: source.number("geom_bond", "distance", "geom_bond_distance", row),
            site_symmetry: source
                .identifier(
                    "geom_bond",
                    "site_symmetry_2",
                    "geom_bond_site_symmetry_2",
                    row,
                )
                .map(Into::into),
            bond_type: source
                .identifier("geom_bond", "type", "geom_bond_type", row)
                .map(Into::into),
        });
    }
    Ok(bonds)
}

fn symmetry_operations(source: &Source<'_>) -> Vec<Box<str>> {
    let column = match source.dialect {
        SmallCifDialect::Ddl2 => source.column(
            "space_group_symop",
            "operation_xyz",
            "space_group_symop_operation_xyz",
        ),
        SmallCifDialect::Ddl1 => {
            source.column("symmetry_equiv", "pos_as_xyz", "symmetry_equiv_pos_as_xyz")
        }
    };
    column
        .into_iter()
        .flat_map(crate::Column::iter)
        .filter_map(CifValue::as_str)
        .map(Into::into)
        .collect()
}
struct Source<'a> {
    block: &'a DataBlock,
    dialect: SmallCifDialect,
}

impl<'a> Source<'a> {
    fn column(&self, category: &str, item: &str, ddl1_tag: &str) -> Option<&'a crate::Column> {
        match self.dialect {
            SmallCifDialect::Ddl2 => self.block.category(category)?.column(item),
            SmallCifDialect::Ddl1 => self.block.category(ddl1_tag)?.column(ddl1_tag),
        }
    }

    fn number(&self, category: &str, item: &str, ddl1_tag: &str, row: usize) -> Option<f64> {
        self.column(category, item, ddl1_tag)?.get(row)?.as_float()
    }

    fn integer(&self, category: &str, item: &str, ddl1_tag: &str, row: usize) -> Option<i32> {
        i32::try_from(
            self.column(category, item, ddl1_tag)?
                .get(row)?
                .as_integer()?,
        )
        .ok()
    }

    fn identifier(
        &self,
        category: &str,
        item: &str,
        ddl1_tag: &str,
        row: usize,
    ) -> Option<std::borrow::Cow<'a, str>> {
        self.column(category, item, ddl1_tag)?
            .get(row)?
            .as_identifier()
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
