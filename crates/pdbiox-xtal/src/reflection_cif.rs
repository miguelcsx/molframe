//! Structure-factor mmCIF lowering and writing over the shared CIF document.

use crate::{
    ReflectionColumn, ReflectionColumnType, ReflectionDataset, ReflectionError, ReflectionTable,
    ReflectionValue,
};
use pdbiox_cif::lexer::Quoting;
use pdbiox_cif::{CifValue, DataBlock, Document, write_preserving};
use pdbiox_core::ByteSpan;
use pdbiox_core::structure::UnitCell;

/// Lowers a lossless mmCIF document's `refln` category to a reflection table.
/// Every `refln` item is retained in source order, including extension items.
///
/// # Errors
///
/// Returns an error for a missing/ragged reflection category or malformed
/// recorded crystal metadata.
pub fn lower_structure_factor_cif(document: &Document) -> Result<ReflectionTable, ReflectionError> {
    let block = document
        .first_block()
        .ok_or(ReflectionError::MissingTable)?;
    let category = block
        .category("refln")
        .ok_or(ReflectionError::MissingTable)?;
    let rows = category.row_count();
    if rows == 0 {
        return Err(ReflectionError::MissingTable);
    }
    let mut columns = Vec::with_capacity(category.len());
    for item in category.items() {
        let source = category.column(item).ok_or(ReflectionError::ColumnLength)?;
        if source.len() != rows {
            return Err(ReflectionError::ColumnLength);
        }
        columns.push(ReflectionColumn {
            label: item.into(),
            column_type: cif_column_type(item),
            values: source.iter().map(from_cif_value).collect(),
            dataset_id: 0,
            mtz_type: None,
        });
    }
    let cell = lower_cell(block)?;
    let space_group_number = integer_item(block, "symmetry", "Int_Tables_number")
        .or_else(|| integer_item(block, "space_group", "IT_number"));
    let space_group_name = text_item(block, "symmetry", "space_group_name_H-M")
        .or_else(|| text_item(block, "space_group", "name_H-M_alt"))
        .map(Into::into);
    let wavelength = float_item(block, "diffrn_radiation_wavelength", "wavelength");
    let datasets = match wavelength {
        Some(value) => vec![ReflectionDataset {
            id: 0,
            project: "".into(),
            crystal: block.name().into(),
            name: block.name().into(),
            wavelength: Some(value),
            cell,
        }],
        None => Vec::new(),
    };
    let table = ReflectionTable {
        title: block.name().into(),
        cell,
        space_group_number,
        space_group_name,
        columns,
        datasets,
        history: Vec::new(),
        symmetry_operations: lower_symmetry_operations(block),
        sort_order: [0; 5],
        resolution_range: None,
        missing_value: None,
        extra_header_records: Vec::new(),
    };
    table.validate()?;
    Ok(table)
}

/// Writes a reflection table as canonical structure-factor mmCIF using the
/// shared lossless CIF writer.
///
/// # Errors
///
/// Returns an error for an invalid/ragged table.
pub fn write_structure_factor_cif(table: &ReflectionTable) -> Result<String, ReflectionError> {
    table.validate()?;
    let name = if table.title.is_empty() {
        "pdbiox"
    } else {
        &table.title
    };
    let mut block = DataBlock::new(name);
    if let Some(cell) = table.cell {
        let category = block.category_mut("cell", ByteSpan::default());
        for (item, value) in [
            ("length_a", cell.lengths[0]),
            ("length_b", cell.lengths[1]),
            ("length_c", cell.lengths[2]),
            ("angle_alpha", cell.angles[0]),
            ("angle_beta", cell.angles[1]),
            ("angle_gamma", cell.angles[2]),
        ] {
            category
                .column_mut(item)
                .push(CifValue::Float(value), Quoting::Bare);
        }
    }
    if table.space_group_number.is_some() || table.space_group_name.is_some() {
        let category = block.category_mut("symmetry", ByteSpan::default());
        if let Some(number) = table.space_group_number {
            category
                .column_mut("Int_Tables_number")
                .push(CifValue::Integer(i64::from(number)), Quoting::Bare);
        }
        if let Some(name) = &table.space_group_name {
            category
                .column_mut("space_group_name_H-M")
                .push(CifValue::Text(name.as_ref().into()), Quoting::Bare);
        }
    }
    if !table.symmetry_operations.is_empty() {
        let category = block.category_mut("space_group_symop", ByteSpan::default());
        for (index, operation) in table.symmetry_operations.iter().enumerate() {
            category.column_mut("id").push(
                CifValue::Integer(
                    i64::try_from(index + 1).map_err(|_| ReflectionError::InvalidNumber)?,
                ),
                Quoting::Bare,
            );
            category
                .column_mut("operation_xyz")
                .push(CifValue::Text(operation.as_ref().into()), Quoting::Bare);
        }
    }
    if let Some(wavelength) = table.datasets.iter().find_map(|dataset| dataset.wavelength) {
        if !wavelength.is_finite() || wavelength <= 0.0 {
            return Err(ReflectionError::InvalidNumber);
        }
        block
            .category_mut("diffrn_radiation_wavelength", ByteSpan::default())
            .column_mut("wavelength")
            .push(CifValue::Float(wavelength), Quoting::Bare);
    }
    let refln = block.category_mut("refln", ByteSpan::default());
    for column in &table.columns {
        let destination = refln.column_mut(&column.label);
        for value in &column.values {
            destination.push(to_cif_value(value), Quoting::Bare);
        }
    }
    let mut document = Document::new();
    document.push(block);
    Ok(write_preserving(&document))
}

fn lower_cell(block: &DataBlock) -> Result<Option<UnitCell>, ReflectionError> {
    let values = [
        float_item(block, "cell", "length_a"),
        float_item(block, "cell", "length_b"),
        float_item(block, "cell", "length_c"),
        float_item(block, "cell", "angle_alpha"),
        float_item(block, "cell", "angle_beta"),
        float_item(block, "cell", "angle_gamma"),
    ];
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
        return Err(ReflectionError::Metadata);
    };
    if [a, b, c]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        || [alpha, beta, gamma]
            .iter()
            .any(|value| !value.is_finite() || !(0.0..180.0).contains(value))
    {
        return Err(ReflectionError::Metadata);
    }
    Ok(Some(UnitCell {
        lengths: [a, b, c],
        angles: [alpha, beta, gamma],
    }))
}

fn float_item(block: &DataBlock, category: &str, item: &str) -> Option<f64> {
    block.category(category)?.value(item, 0)?.as_float()
}
fn integer_item(block: &DataBlock, category: &str, item: &str) -> Option<i32> {
    i32::try_from(block.category(category)?.value(item, 0)?.as_integer()?).ok()
}
fn text_item<'a>(block: &'a DataBlock, category: &str, item: &str) -> Option<&'a str> {
    block.category(category)?.text(item, 0)
}

fn lower_symmetry_operations(block: &DataBlock) -> Vec<Box<str>> {
    let Some(category) = block
        .category("space_group_symop")
        .or_else(|| block.category("symmetry_equiv"))
    else {
        return Vec::new();
    };
    let item = if category.column("operation_xyz").is_some() {
        "operation_xyz"
    } else {
        "pos_as_xyz"
    };
    category
        .column(item)
        .into_iter()
        .flat_map(pdbiox_cif::Column::iter)
        .filter_map(CifValue::as_str)
        .map(Into::into)
        .collect()
}

fn from_cif_value(value: &CifValue) -> ReflectionValue {
    match value {
        CifValue::Unknown => ReflectionValue::Missing,
        CifValue::Inapplicable => ReflectionValue::Inapplicable,
        CifValue::Integer(value) => ReflectionValue::Integer(*value),
        CifValue::Float(value) => ReflectionValue::Real(*value),
        CifValue::Text(value) => ReflectionValue::Text(value.as_ref().into()),
    }
}
fn to_cif_value(value: &ReflectionValue) -> CifValue {
    match value {
        ReflectionValue::Missing => CifValue::Unknown,
        ReflectionValue::Inapplicable => CifValue::Inapplicable,
        ReflectionValue::Integer(value) => CifValue::Integer(*value),
        ReflectionValue::Real(value) => CifValue::Float(*value),
        ReflectionValue::Text(value) => CifValue::Text(value.as_ref().into()),
    }
}

fn cif_column_type(item: &str) -> ReflectionColumnType {
    let lower = item.to_ascii_lowercase();
    if matches!(lower.as_str(), "index_h" | "index_k" | "index_l") {
        ReflectionColumnType::MillerIndex
    } else if lower.contains("sigma") {
        ReflectionColumnType::StandardDeviation
    } else if lower.contains("phase") {
        ReflectionColumnType::Phase
    } else if lower.contains("squared") || lower.contains("intensity") {
        ReflectionColumnType::Intensity
    } else if lower.contains("f_meas") || lower.contains("f_calc") {
        ReflectionColumnType::Amplitude
    } else if lower.contains("status") || lower.contains("details") {
        ReflectionColumnType::Text
    } else if lower.contains("free") || lower.ends_with("flag") {
        ReflectionColumnType::Flag
    } else {
        ReflectionColumnType::Real
    }
}

#[cfg(test)]
#[path = "reflection_cif_tests.rs"]
mod tests;
