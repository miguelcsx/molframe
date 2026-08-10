//! Merged CCP4 MTZ reflection reader and writer.

use crate::{
    CellTransform, ReflectionColumn, ReflectionColumnType, ReflectionDataset, ReflectionError,
    ReflectionTable, ReflectionValue,
};
use pdbiox_core::structure::UnitCell;
use std::collections::BTreeMap;

use crate::numeric::{f32_to_i64, f64_to_f32, i64_to_f32};

#[path = "mtz_header_parse.rs"]
mod header_parse;
use header_parse::{parse_header, parse_history};
#[path = "mtz_header_write.rs"]
mod header_write;
use header_write::{
    write_column_headers, write_crystal_header, write_dataset_headers, write_header_tail,
};

const PREAMBLE_BYTES: usize = 80;
const RECORD_BYTES: usize = 80;

#[derive(Clone, Copy)]
enum Endian {
    Little,
    Big,
}

#[derive(Default)]
struct Header {
    title: Box<str>,
    columns: Vec<ColumnHeader>,
    reflection_count: usize,
    batch_count: usize,
    cell: Option<UnitCell>,
    sort_order: [i32; 5],
    space_group_number: Option<i32>,
    space_group_name: Option<Box<str>>,
    symmetry_operations: Vec<Box<str>>,
    resolution_range: Option<[f64; 2]>,
    missing_value: Option<f32>,
    datasets: BTreeMap<i32, DatasetHeader>,
    extras: Vec<Box<str>>,
    next_offset: usize,
}

struct ColumnHeader {
    label: Box<str>,
    kind: char,
    dataset_id: i32,
}

#[derive(Default)]
struct DatasetHeader {
    project: Box<str>,
    crystal: Box<str>,
    name: Box<str>,
    wavelength: Option<f64>,
    cell: Option<UnitCell>,
}

/// Reads a merged MTZ file, preserving columns, datasets, symmetry, missing
/// values, history, and unknown main-header records.
///
/// Unmerged batch MTZ files are refused explicitly because interpreting their
/// batch binary headers as merged observations would corrupt semantics.
///
/// # Errors
///
/// Returns an error for malformed/truncated framing, inconsistent counts,
/// unsupported batch data, or invalid numbers.
pub fn read_mtz(bytes: &[u8]) -> Result<ReflectionTable, ReflectionError> {
    if bytes.len() < PREAMBLE_BYTES || &bytes[..4] != b"MTZ " {
        return Err(ReflectionError::TruncatedMtz);
    }
    let endian = detect_endian(bytes)?;
    let header_word = read_header_word(bytes, endian)?;
    let header_offset = header_word
        .checked_sub(1)
        .and_then(|word| usize::try_from(word).ok())
        .and_then(|word| word.checked_mul(4))
        .ok_or(ReflectionError::InvalidMtz)?;
    if header_offset < PREAMBLE_BYTES || header_offset > bytes.len() {
        return Err(ReflectionError::InvalidMtz);
    }
    let header = parse_header(bytes, header_offset)?;
    if header.batch_count != 0 {
        return Err(ReflectionError::Unsupported(
            "unmerged MTZ batch headers".into(),
        ));
    }
    if header.columns.is_empty() || header.reflection_count == 0 {
        return Err(ReflectionError::MissingTable);
    }
    let value_count = header
        .columns
        .len()
        .checked_mul(header.reflection_count)
        .ok_or(ReflectionError::InvalidMtz)?;
    let expected_header = PREAMBLE_BYTES
        .checked_add(
            value_count
                .checked_mul(4)
                .ok_or(ReflectionError::InvalidMtz)?,
        )
        .ok_or(ReflectionError::InvalidMtz)?;
    if expected_header != header_offset {
        return Err(ReflectionError::InvalidMtz);
    }
    let mut columns: Vec<_> = header
        .columns
        .iter()
        .map(|column| ReflectionColumn {
            label: column.label.clone(),
            column_type: mtz_column_type(column.kind),
            values: Vec::with_capacity(header.reflection_count),
            dataset_id: column.dataset_id,
            mtz_type: Some(column.kind),
        })
        .collect();
    let column_count = columns.len();
    for row in 0..header.reflection_count {
        for (column_index, column) in columns.iter_mut().enumerate() {
            let index = row * column_count + column_index;
            let value = read_f32(bytes, PREAMBLE_BYTES + index * 4, endian)?;
            let missing = value.is_nan()
                || header
                    .missing_value
                    .is_some_and(|sentinel| value.to_bits() == sentinel.to_bits());
            column.values.push(if missing {
                ReflectionValue::Missing
            } else if matches!(
                column.column_type,
                ReflectionColumnType::MillerIndex | ReflectionColumnType::Flag
            ) && value.is_finite()
                && value.fract().abs() <= f32::EPSILON
            {
                ReflectionValue::Integer(f32_to_i64(value))
            } else {
                ReflectionValue::Real(f64::from(value))
            });
        }
    }
    let (history, _) = parse_history(bytes, header.next_offset)?;
    let datasets = header
        .datasets
        .into_iter()
        .map(|(id, dataset)| ReflectionDataset {
            id,
            project: dataset.project,
            crystal: dataset.crystal,
            name: dataset.name,
            wavelength: dataset.wavelength,
            cell: dataset.cell,
        })
        .collect();
    let table = ReflectionTable {
        title: header.title,
        cell: header.cell,
        space_group_number: header.space_group_number,
        space_group_name: header.space_group_name,
        columns,
        datasets,
        history,
        symmetry_operations: header.symmetry_operations,
        sort_order: header.sort_order,
        resolution_range: header.resolution_range,
        missing_value: header.missing_value,
        extra_header_records: header.extras,
    };
    table.validate()?;
    Ok(table)
}

/// Writes a merged, little-endian MTZ V1.1 file.
///
/// # Errors
///
/// Returns an error for invalid tables, absent required crystal/symmetry
/// metadata, text columns, labels beyond MTZ limits, or values not exactly
/// representable as MTZ `REAL*4` data.
pub fn write_mtz(table: &ReflectionTable) -> Result<Vec<u8>, ReflectionError> {
    table.validate()?;
    let cell = table.cell.ok_or(ReflectionError::Metadata)?;
    let space_group_number = table.space_group_number.ok_or(ReflectionError::Metadata)?;
    let space_group_name = table
        .space_group_name
        .as_deref()
        .ok_or(ReflectionError::Metadata)?;
    if table.symmetry_operations.is_empty() {
        return Err(ReflectionError::Unsupported(
            "MTZ requires explicit symmetry operations".into(),
        ));
    }
    if table.columns.len() < 3 || table.columns.iter().any(|column| column.label.len() > 30) {
        return Err(ReflectionError::InvalidMtz);
    }
    let value_count = table
        .row_count()
        .checked_mul(table.columns.len())
        .ok_or(ReflectionError::InvalidMtz)?;
    let header_word = value_count
        .checked_add(21)
        .and_then(|value| i32::try_from(value).ok())
        .ok_or(ReflectionError::Unsupported(
            "64-bit MTZ header offsets".into(),
        ))?;
    let mut output = vec![0_u8; PREAMBLE_BYTES];
    output[..4].copy_from_slice(b"MTZ ");
    output[4..8].copy_from_slice(&header_word.to_le_bytes());
    output[8..12].copy_from_slice(&[0x44, 0x41, 0, 0]);
    for row in 0..table.row_count() {
        for column in &table.columns {
            let value = encode_value(&column.values[row], table.missing_value)?;
            output.extend_from_slice(&value.to_le_bytes());
        }
    }
    write_crystal_header(
        &mut output,
        table,
        cell,
        space_group_number,
        space_group_name,
    )?;
    write_column_headers(&mut output, table)?;
    write_dataset_headers(&mut output, table, cell)?;
    write_header_tail(&mut output, table);
    Ok(output)
}

fn detect_endian(bytes: &[u8]) -> Result<Endian, ReflectionError> {
    match &bytes[8..10] {
        [0x44, 0x41 | 0x44] => Ok(Endian::Little),
        [0x11, 0x11] => Ok(Endian::Big),
        _ => Err(ReflectionError::InvalidMtz),
    }
}
fn read_header_word(bytes: &[u8], endian: Endian) -> Result<i64, ReflectionError> {
    let value = read_i32(bytes, 4, endian)?;
    if value != -1 {
        return Ok(i64::from(value));
    }
    let raw: [u8; 8] = bytes
        .get(12..20)
        .ok_or(ReflectionError::TruncatedMtz)?
        .try_into()
        .map_err(|_| ReflectionError::TruncatedMtz)?;
    Ok(match endian {
        Endian::Little => i64::from_le_bytes(raw),
        Endian::Big => i64::from_be_bytes(raw),
    })
}
fn read_i32(bytes: &[u8], offset: usize, endian: Endian) -> Result<i32, ReflectionError> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or(ReflectionError::TruncatedMtz)?
        .try_into()
        .map_err(|_| ReflectionError::TruncatedMtz)?;
    Ok(match endian {
        Endian::Little => i32::from_le_bytes(raw),
        Endian::Big => i32::from_be_bytes(raw),
    })
}
fn read_f32(bytes: &[u8], offset: usize, endian: Endian) -> Result<f32, ReflectionError> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or(ReflectionError::TruncatedMtz)?
        .try_into()
        .map_err(|_| ReflectionError::TruncatedMtz)?;
    Ok(match endian {
        Endian::Little => f32::from_le_bytes(raw),
        Endian::Big => f32::from_be_bytes(raw),
    })
}
fn record(bytes: &[u8], offset: usize) -> Result<&str, ReflectionError> {
    let raw = bytes
        .get(offset..offset + RECORD_BYTES)
        .ok_or(ReflectionError::TruncatedMtz)?;
    std::str::from_utf8(raw).map_err(|_| ReflectionError::InvalidMtz)
}
fn split_keyword(line: &str) -> (&str, &str) {
    line.split_once(char::is_whitespace)
        .map_or((line, ""), |(key, rest)| (key, rest.trim_start()))
}
fn next_parse<'a, T: std::str::FromStr>(
    fields: &mut impl Iterator<Item = &'a str>,
) -> Result<T, ReflectionError> {
    fields
        .next()
        .ok_or(ReflectionError::InvalidMtz)?
        .parse()
        .map_err(|_| ReflectionError::InvalidMtz)
}
fn parse_i64_array<const N: usize>(text: &str) -> Result<[i64; N], ReflectionError> {
    let values = text
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| ReflectionError::InvalidMtz)?;
    values.try_into().map_err(|_| ReflectionError::InvalidMtz)
}
fn parse_i32_array<const N: usize>(text: &str) -> Result<[i32; N], ReflectionError> {
    let values = text
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| ReflectionError::InvalidMtz)?;
    values.try_into().map_err(|_| ReflectionError::InvalidMtz)
}
fn parse_f64_array<const N: usize>(text: &str) -> Result<[f64; N], ReflectionError> {
    let values = text
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| ReflectionError::InvalidMtz)?;
    values.try_into().map_err(|_| ReflectionError::InvalidMtz)
}
fn parse_cell(text: &str) -> Result<UnitCell, ReflectionError> {
    let values = parse_f64_array::<6>(text)?;
    if values[..3]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        || values[3..]
            .iter()
            .any(|value| !value.is_finite() || !(0.0..180.0).contains(value))
    {
        return Err(ReflectionError::Metadata);
    }
    Ok(UnitCell {
        lengths: [values[0], values[1], values[2]],
        angles: [values[3], values[4], values[5]],
    })
}
fn parse_missing(text: &str) -> Result<Option<f32>, ReflectionError> {
    if text.eq_ignore_ascii_case("NAN") {
        Ok(None)
    } else {
        let value = text.parse().map_err(|_| ReflectionError::InvalidMtz)?;
        if f32::is_finite(value) {
            Ok(Some(value))
        } else {
            Err(ReflectionError::InvalidNumber)
        }
    }
}
fn mtz_column_type(kind: char) -> ReflectionColumnType {
    match kind {
        'H' => ReflectionColumnType::MillerIndex,
        'F' | 'D' | 'G' | 'E' => ReflectionColumnType::Amplitude,
        'J' | 'K' => ReflectionColumnType::Intensity,
        'Q' | 'L' | 'M' => ReflectionColumnType::StandardDeviation,
        'P' => ReflectionColumnType::Phase,
        'B' | 'Y' | 'I' => ReflectionColumnType::Flag,
        _ => ReflectionColumnType::Real,
    }
}
fn mtz_type(kind: ReflectionColumnType) -> char {
    match kind {
        ReflectionColumnType::MillerIndex => 'H',
        ReflectionColumnType::Amplitude => 'F',
        ReflectionColumnType::Intensity => 'J',
        ReflectionColumnType::StandardDeviation => 'Q',
        ReflectionColumnType::Phase => 'P',
        ReflectionColumnType::Flag => 'I',
        ReflectionColumnType::Real => 'R',
        ReflectionColumnType::Text => 'S',
    }
}
fn encode_value(value: &ReflectionValue, missing: Option<f32>) -> Result<f32, ReflectionError> {
    match value {
        ReflectionValue::Missing | ReflectionValue::Inapplicable => match missing {
            Some(value) => Ok(value),
            None => Ok(f32::NAN),
        },
        ReflectionValue::Integer(value) => {
            let encoded = i64_to_f32(*value);
            if encoded.is_finite() && f32_to_i64(encoded) == *value {
                Ok(encoded)
            } else {
                Err(ReflectionError::InvalidNumber)
            }
        }
        ReflectionValue::Real(value) => {
            let encoded = f64_to_f32(*value);
            if value.is_finite() && encoded.is_finite() {
                Ok(encoded)
            } else {
                Err(ReflectionError::InvalidNumber)
            }
        }
        ReflectionValue::Text(_) => Err(ReflectionError::Unsupported(
            "text-valued MTZ column".into(),
        )),
    }
}
fn column_min_max(column: &ReflectionColumn) -> Result<(f32, f32), ReflectionError> {
    let mut minimum = f32::INFINITY;
    let mut maximum = f32::NEG_INFINITY;
    for value in &column.values {
        if matches!(
            value,
            ReflectionValue::Missing | ReflectionValue::Inapplicable
        ) {
            continue;
        }
        let encoded = encode_value(value, None)?;
        minimum = minimum.min(encoded);
        maximum = maximum.max(encoded);
    }
    if minimum.is_infinite() {
        Ok((0.0, 0.0))
    } else {
        Ok((minimum, maximum))
    }
}
fn format_cell(keyword: &str, id: Option<i32>, cell: UnitCell) -> String {
    let prefix = id.map_or_else(|| keyword.to_owned(), |id| format!("{keyword} {id}"));
    format!(
        "{prefix} {:.5} {:.5} {:.5} {:.5} {:.5} {:.5}",
        cell.lengths[0],
        cell.lengths[1],
        cell.lengths[2],
        cell.angles[0],
        cell.angles[1],
        cell.angles[2]
    )
}
fn push_record(output: &mut Vec<u8>, text: &str) {
    let source = text.as_bytes();
    let length = source.len().min(RECORD_BYTES);
    output.extend_from_slice(&source[..length]);
    output.resize(output.len() + RECORD_BYTES - length, b' ');
}
fn calculate_resolution(
    table: &ReflectionTable,
    cell: UnitCell,
) -> Result<[f64; 2], ReflectionError> {
    let transform = CellTransform::new(&cell).map_err(|_| ReflectionError::Metadata)?;
    let inverse = transform.inverse_matrix();
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for hkl in table.miller_indices()? {
        let reciprocal = (0..3)
            .map(|axis| {
                (0..3)
                    .map(|index| inverse[index][axis] * f64::from(hkl[index]))
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let value = reciprocal
            .iter()
            .map(|component| component * component)
            .sum::<f64>();
        minimum = minimum.min(value);
        maximum = maximum.max(value);
    }
    Ok([minimum, maximum])
}

#[cfg(test)]
#[path = "mtz_tests.rs"]
mod tests;
