//! Central DMS schema vocabulary and structural validation.

use crate::dms::DmsError;
use rusqlite::Connection;
use std::collections::{BTreeMap, BTreeSet};

pub(crate) const PARTICLE: &str = "particle";
pub(crate) const BOND: &str = "bond";
pub(crate) const GLOBAL_CELL: &str = "global_cell";
pub(crate) const VERSION: &str = "dms_version";

pub(crate) const REQUIRED_PARTICLE: [&str; 4] = ["id", "x", "y", "z"];
pub(crate) const REQUIRED_BOND: [&str; 3] = ["p0", "p1", "order"];
pub(crate) const REQUIRED_CELL: [&str; 4] = ["id", "x", "y", "z"];

pub(crate) const PARTICLE_PROJECTION: [&str; 25] = [
    "id",
    "anum",
    "msys_ct",
    "nbtype",
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "mass",
    "charge",
    "resid",
    "resname",
    "chain",
    "segid",
    "name",
    "insertion",
    "formal_charge",
    "occupancy",
    "bfactor",
    "grp_temperature",
    "grp_energy",
    "grp_ligand",
    "grp_bias",
];

pub(crate) const PARTICLE_INTEGERS: [&str; 9] = [
    "id",
    "anum",
    "msys_ct",
    "nbtype",
    "resid",
    "grp_temperature",
    "grp_energy",
    "grp_ligand",
    "grp_bias",
];
pub(crate) const PARTICLE_FLOATS: [&str; 11] = [
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "mass",
    "charge",
    "formal_charge",
    "occupancy",
    "bfactor",
];
pub(crate) const PARTICLE_TEXT: [&str; 5] = ["resname", "chain", "segid", "name", "insertion"];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Affinity {
    Integer,
    Float,
    Text,
}

pub(crate) fn table_columns(
    connection: &Connection,
    table: &'static str,
) -> Result<BTreeSet<Box<str>>, DmsError> {
    if !table_exists(connection, table)? {
        return Err(DmsError::MissingTable { table });
    }
    let sql = format!("PRAGMA table_info({table})");
    let mut statement = connection.prepare(&sql)?;
    let rows = statement.query_map([], |row| row.get::<_, String>(1))?;
    let mut columns = BTreeSet::new();
    for row in rows {
        columns.insert(row?.to_ascii_lowercase().into_boxed_str());
    }
    Ok(columns)
}

pub(crate) fn optional_table_columns(
    connection: &Connection,
    table: &'static str,
) -> Result<Option<BTreeSet<Box<str>>>, DmsError> {
    if table_exists(connection, table)? {
        table_columns(connection, table).map(Some)
    } else {
        Ok(None)
    }
}

pub(crate) fn require_columns(
    table: &'static str,
    columns: &BTreeSet<Box<str>>,
    required: &[&'static str],
) -> Result<(), DmsError> {
    for column in required {
        if !columns.contains(*column) {
            return Err(DmsError::MissingColumn { table, column });
        }
    }
    Ok(())
}

pub(crate) fn validate_affinities(
    connection: &Connection,
    table: &'static str,
    expected: &[(&'static str, Affinity)],
) -> Result<(), DmsError> {
    let sql = format!("PRAGMA table_info({table})");
    let mut statement = connection.prepare(&sql)?;
    let rows = statement.query_map([], |row| {
        Ok((
            row.get::<_, String>(1)?.to_ascii_lowercase(),
            row.get::<_, String>(2)?,
        ))
    })?;
    let declared = rows.collect::<Result<BTreeMap<_, _>, _>>()?;
    for (column, expected_affinity) in expected {
        let Some(declared_type) = declared.get(*column) else {
            continue;
        };
        if !declared_type.trim().is_empty() && affinity(declared_type) != Some(*expected_affinity) {
            return Err(DmsError::InvalidColumnType { table, column });
        }
    }
    Ok(())
}

pub(crate) fn validate_standard_types(
    connection: &Connection,
    table: &'static str,
    integers: &[&'static str],
    floats: &[&'static str],
    text: &[&'static str],
) -> Result<(), DmsError> {
    let expected = integers
        .iter()
        .map(|column| (*column, Affinity::Integer))
        .chain(floats.iter().map(|column| (*column, Affinity::Float)))
        .chain(text.iter().map(|column| (*column, Affinity::Text)))
        .collect::<Vec<_>>();
    validate_affinities(connection, table, &expected)
}

pub(crate) fn select_projection(
    columns: &BTreeSet<Box<str>>,
    projection: &[&'static str],
) -> String {
    projection
        .iter()
        .map(|column| {
            if columns.contains(*column) {
                (*column).to_owned()
            } else {
                format!("NULL AS {column}")
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn table_exists(connection: &Connection, table: &'static str) -> Result<bool, DmsError> {
    let mut statement = connection.prepare(
        "SELECT 1 FROM sqlite_schema WHERE type IN ('table', 'view') AND lower(name) = ?1 LIMIT 1",
    )?;
    Ok(statement.exists([table])?)
}

fn affinity(declared: &str) -> Option<Affinity> {
    let declared = declared.to_ascii_uppercase();
    if declared.contains("INT") {
        Some(Affinity::Integer)
    } else if declared.contains("CHAR") || declared.contains("CLOB") || declared.contains("TEXT") {
        Some(Affinity::Text)
    } else if declared.contains("REAL") || declared.contains("FLOA") || declared.contains("DOUB") {
        Some(Affinity::Float)
    } else {
        None
    }
}
