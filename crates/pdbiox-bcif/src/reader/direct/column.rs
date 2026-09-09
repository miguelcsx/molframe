//! Compact decoded coordinate columns and borrowed row access.

use super::container::{DirectColumn, DirectData};
use super::values::{ColumnValues, ValueRef};
use crate::codec::{Decoded, decode_borrowed, decode_f32_borrowed_into};
use num_traits::ToPrimitive;
use pdbiox_cif::{AtomSiteRow, AtomSiteRowSink, Field};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use std::borrow::Cow;

pub(super) const FIELD_COUNT: usize = 21;

pub(super) struct AtomColumns {
    columns: [Option<DecodedColumn>; FIELD_COUNT],
    row_count: usize,
}

impl AtomColumns {
    pub(super) fn new(row_count: usize) -> Self {
        Self {
            columns: std::array::from_fn(|_| None),
            row_count,
        }
    }

    pub(super) fn insert(
        &mut self,
        item: &str,
        encoded: DirectColumn<'_>,
    ) -> Result<(), Diagnostic> {
        let Some(index) = field_index(item) else {
            return Ok(());
        };
        if self.columns[index].is_some() {
            return Err(duplicate_item(item));
        }
        self.columns[index] = Some(DecodedColumn::decode(
            encoded,
            self.row_count,
            (14..=18).contains(&index),
        )?);
        Ok(())
    }

    pub(super) fn column(&self, item: &str) -> Option<&DecodedColumn> {
        self.columns.get(field_index(item)?)?.as_ref()
    }

    pub(super) fn row_count(&self) -> usize {
        self.row_count
    }

    pub(super) fn feed(&self, sink: &mut dyn AtomSiteRowSink) {
        for row in 0..self.row_count {
            sink.feed(&BinaryAtomRow { columns: self, row });
        }
    }

    pub(super) fn compact_floats_for_lowering(&mut self) {
        for index in 14..=18 {
            if let Some(column) = &mut self.columns[index] {
                column.values.compact_float();
            }
        }
    }

    pub(super) fn iter(&self) -> impl Iterator<Item = (&'static str, &DecodedColumn)> {
        FIELD_NAMES
            .iter()
            .copied()
            .zip(&self.columns)
            .filter_map(|(name, column)| column.as_ref().map(|column| (name, column)))
    }
}

pub(super) struct DecodedColumn {
    values: ColumnValues,
    mask: Option<Vec<u8>>,
}

impl DecodedColumn {
    pub(super) fn decode(
        encoded: DirectColumn<'_>,
        row_count: usize,
        compact_float: bool,
    ) -> Result<Self, Diagnostic> {
        let values = if compact_float {
            let mut values = vec![0.0_f32; row_count];
            decode_f32_borrowed_into(&encoded.data.encoding, encoded.data.data, &mut values)?;
            ColumnValues::from_f32(values)
        } else {
            ColumnValues::new(decode_borrowed(&encoded.data.encoding, encoded.data.data)?)
        };
        if values.len() != row_count {
            return Err(length_error(values.len(), row_count));
        }
        let mask = match encoded.mask {
            Some(encoded) => Some(decode_mask(&encoded, row_count)?),
            None => None,
        };
        Ok(Self { values, mask })
    }

    pub(super) fn integer(&self, row: usize) -> Option<i64> {
        match self.value(row)? {
            ValueRef::Integer(value) => Some(value),
            _ => None,
        }
    }

    pub(super) fn same_rows(&self, left: usize, right: usize) -> bool {
        match (self.value(left), self.value(right)) {
            (Some(ValueRef::Inapplicable), Some(ValueRef::Inapplicable))
            | (Some(ValueRef::Unknown), Some(ValueRef::Unknown)) => true,
            (Some(ValueRef::Text(left)), Some(ValueRef::Text(right))) => left == right,
            (Some(ValueRef::Integer(left)), Some(ValueRef::Integer(right))) => left == right,
            (Some(ValueRef::Float(left)), Some(ValueRef::Float(right))) => {
                left.partial_cmp(&right) == Some(std::cmp::Ordering::Equal)
            }
            _ => false,
        }
    }

    fn value(&self, row: usize) -> Option<ValueRef<'_>> {
        match self.mask.as_deref().and_then(|mask| mask.get(row)).copied() {
            Some(1) => return Some(ValueRef::Inapplicable),
            Some(2) => return Some(ValueRef::Unknown),
            Some(0) | None => {}
            Some(_) => return None,
        }
        self.values.value(row)
    }
}

struct BinaryAtomRow<'a> {
    columns: &'a AtomColumns,
    row: usize,
}

impl BinaryAtomRow<'_> {
    fn value(&self, field: Field) -> Option<ValueRef<'_>> {
        self.columns.column(field.item())?.value(self.row)
    }
}

impl AtomSiteRow for BinaryAtomRow<'_> {
    fn row(&self) -> usize {
        self.row
    }

    fn text(&self, field: Field) -> Option<&str> {
        match self.value(field)? {
            ValueRef::Text(value) => Some(value),
            _ => None,
        }
    }

    fn identifier(&self, field: Field) -> Option<Cow<'_, str>> {
        match self.value(field)? {
            ValueRef::Text(value) => Some(Cow::Borrowed(value)),
            ValueRef::Integer(value) => Some(Cow::Owned(value.to_string())),
            ValueRef::Float(value) => Some(Cow::Owned(value.to_string())),
            ValueRef::Inapplicable | ValueRef::Unknown => None,
        }
    }

    fn integer(&self, field: Field) -> Option<i64> {
        match self.value(field)? {
            ValueRef::Integer(value) => Some(value),
            _ => None,
        }
    }

    fn float(&self, field: Field) -> Option<f64> {
        match self.value(field)? {
            ValueRef::Integer(value) => value.to_f64(),
            ValueRef::Float(value) => Some(value),
            _ => None,
        }
    }

    fn is_recorded(&self, field: Field) -> bool {
        matches!(
            self.value(field),
            Some(ValueRef::Text(_) | ValueRef::Integer(_) | ValueRef::Float(_))
        )
    }
}

pub(super) fn field_index(item: &str) -> Option<usize> {
    match item.as_bytes() {
        b"group_PDB" => Some(0),
        b"id" => Some(1),
        b"type_symbol" => Some(2),
        b"label_atom_id" => Some(3),
        b"auth_atom_id" => Some(4),
        b"label_alt_id" => Some(5),
        b"label_comp_id" => Some(6),
        b"auth_comp_id" => Some(7),
        b"label_asym_id" => Some(8),
        b"auth_asym_id" => Some(9),
        b"label_entity_id" => Some(10),
        b"label_seq_id" => Some(11),
        b"auth_seq_id" => Some(12),
        b"pdbx_PDB_ins_code" => Some(13),
        b"Cartn_x" => Some(14),
        b"Cartn_y" => Some(15),
        b"Cartn_z" => Some(16),
        b"occupancy" => Some(17),
        b"B_iso_or_equiv" => Some(18),
        b"pdbx_formal_charge" => Some(19),
        b"pdbx_PDB_model_num" => Some(20),
        _ => None,
    }
}

pub(super) fn is_frame_coordinate(item: &str) -> bool {
    matches!(
        item,
        "Cartn_x" | "Cartn_y" | "Cartn_z" | "pdbx_PDB_model_num"
    )
}

fn decode_mask(encoded: &DirectData<'_>, row_count: usize) -> Result<Vec<u8>, Diagnostic> {
    let Decoded::Integers(values) = decode_borrowed(&encoded.encoding, encoded.data)? else {
        return Err(
            Diagnostic::new(Code::E1403).with_message("a BinaryCIF mask is not an integer array")
        );
    };
    if values.len() != row_count {
        return Err(length_error(values.len(), row_count));
    }
    let mut compact = Vec::with_capacity(row_count);
    for value in values {
        let mask = match value {
            0 => 0,
            1 => 1,
            2 => 2,
            _ => {
                return Err(
                    Diagnostic::new(Code::E1401).with_context("mask value", value.to_string())
                );
            }
        };
        compact.push(mask);
    }
    Ok(compact)
}

fn duplicate_item(item: &str) -> Diagnostic {
    Diagnostic::new(Code::E1401)
        .with_message("a BinaryCIF category contains a duplicate item")
        .with_context("item", item)
}

fn length_error(actual: usize, declared: usize) -> Diagnostic {
    Diagnostic::new(Code::E1401)
        .with_context("actual", actual.to_string())
        .with_context("declared", declared.to_string())
}

const FIELD_NAMES: [&str; FIELD_COUNT] = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "auth_atom_id",
    "label_alt_id",
    "label_comp_id",
    "auth_comp_id",
    "label_asym_id",
    "auth_asym_id",
    "label_entity_id",
    "label_seq_id",
    "auth_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "pdbx_formal_charge",
    "pdbx_PDB_model_num",
];
