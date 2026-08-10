//! MTZ textual header parsing.

use super::{
    BTreeMap, ColumnHeader, DatasetHeader, Header, RECORD_BYTES, ReflectionError, next_parse,
    parse_cell, parse_f64_array, parse_i32_array, parse_i64_array, parse_missing, record,
    split_keyword,
};

pub(super) fn parse_header(bytes: &[u8], mut offset: usize) -> Result<Header, ReflectionError> {
    let mut output = Header::default();
    loop {
        let line = record(bytes, offset)?;
        offset += RECORD_BYTES;
        let trimmed = line.trim_end();
        let (keyword, arguments) = split_keyword(trimmed);
        match keyword {
            "END" => break,
            "VERS" | "NDIF" | "" => {}
            "TITLE" => output.title = arguments.into(),
            "NCOL" => parse_ncol(arguments, &mut output)?,
            "CELL" => output.cell = Some(parse_cell(arguments)?),
            "SORT" => output.sort_order = parse_i32_array::<5>(arguments)?,
            "SYMINF" => parse_syminf(arguments, &mut output)?,
            "SYMM" => output.symmetry_operations.push(arguments.into()),
            "RESO" => output.resolution_range = Some(parse_f64_array::<2>(arguments)?),
            "VALM" => output.missing_value = parse_missing(arguments)?,
            "COLUMN" => output.columns.push(parse_column(arguments)?),
            "PROJECT" => update_dataset(arguments, &mut output.datasets, DatasetField::Project)?,
            "CRYSTAL" => update_dataset(arguments, &mut output.datasets, DatasetField::Crystal)?,
            "DATASET" => update_dataset(arguments, &mut output.datasets, DatasetField::Name)?,
            "DCELL" => update_dataset(arguments, &mut output.datasets, DatasetField::Cell)?,
            "DWAVEL" => update_dataset(arguments, &mut output.datasets, DatasetField::Wavelength)?,
            _ => output.extras.push(trimmed.into()),
        }
    }
    if output.columns.len() < 3 || output.columns.len() > 200 {
        return Err(ReflectionError::InvalidMtz);
    }
    output.next_offset = offset;
    Ok(output)
}

fn parse_ncol(arguments: &str, header: &mut Header) -> Result<(), ReflectionError> {
    let values = parse_i64_array::<3>(arguments)?;
    header.reflection_count =
        usize::try_from(values[1]).map_err(|_| ReflectionError::InvalidMtz)?;
    header.batch_count = usize::try_from(values[2]).map_err(|_| ReflectionError::InvalidMtz)?;
    if usize::try_from(values[0])
        .ok()
        .is_some_and(|count| count > 200)
    {
        return Err(ReflectionError::InvalidMtz);
    }
    Ok(())
}

fn parse_syminf(arguments: &str, header: &mut Header) -> Result<(), ReflectionError> {
    let mut fields = arguments.split_whitespace();
    let _operations: usize = next_parse(&mut fields)?;
    let _primitive: usize = next_parse(&mut fields)?;
    let _lattice = fields.next().ok_or(ReflectionError::InvalidMtz)?;
    header.space_group_number = Some(next_parse(&mut fields)?);
    let start = arguments.find('\'').ok_or(ReflectionError::InvalidMtz)? + 1;
    let end = arguments[start..]
        .find('\'')
        .ok_or(ReflectionError::InvalidMtz)?
        + start;
    header.space_group_name = Some(arguments[start..end].into());
    Ok(())
}

fn parse_column(arguments: &str) -> Result<ColumnHeader, ReflectionError> {
    let mut fields = arguments.split_whitespace();
    let label = fields.next().ok_or(ReflectionError::InvalidMtz)?;
    let kind = fields
        .next()
        .and_then(|value| value.chars().next())
        .filter(|_| {
            fields
                .next()
                .and_then(|value| value.parse::<f32>().ok())
                .is_some()
        })
        .ok_or(ReflectionError::InvalidMtz)?;
    let _maximum: f32 = next_parse(&mut fields)?;
    Ok(ColumnHeader {
        label: label.into(),
        kind,
        dataset_id: next_parse(&mut fields)?,
    })
}

#[derive(Clone, Copy)]
enum DatasetField {
    Project,
    Crystal,
    Name,
    Cell,
    Wavelength,
}

fn update_dataset(
    arguments: &str,
    datasets: &mut BTreeMap<i32, DatasetHeader>,
    field: DatasetField,
) -> Result<(), ReflectionError> {
    let mut fields = arguments.split_whitespace();
    let id = next_parse(&mut fields)?;
    let dataset = datasets.entry(id).or_default();
    match field {
        DatasetField::Project => {
            dataset.project = fields.next().ok_or(ReflectionError::InvalidMtz)?.into();
        }
        DatasetField::Crystal => {
            dataset.crystal = fields.next().ok_or(ReflectionError::InvalidMtz)?.into();
        }
        DatasetField::Name => {
            dataset.name = fields.next().ok_or(ReflectionError::InvalidMtz)?.into();
        }
        DatasetField::Cell => {
            dataset.cell = Some(parse_cell(&fields.collect::<Vec<_>>().join(" "))?);
        }
        DatasetField::Wavelength => {
            let value: f64 = next_parse(&mut fields)?;
            if !value.is_finite() {
                return Err(ReflectionError::InvalidNumber);
            }
            dataset.wavelength = Some(value);
        }
    }
    Ok(())
}

pub(super) fn parse_history(
    bytes: &[u8],
    mut offset: usize,
) -> Result<(Vec<Box<str>>, usize), ReflectionError> {
    let mut history = Vec::new();
    while offset + RECORD_BYTES <= bytes.len() {
        let line = record(bytes, offset)?;
        offset += RECORD_BYTES;
        let trimmed = line.trim_end();
        if trimmed.starts_with("MTZENDOFHEADERS") {
            return Ok((history, offset));
        }
        if let Some(count) = trimmed
            .strip_prefix("MTZHIST")
            .and_then(|text| text.trim().parse::<usize>().ok())
        {
            if count > 30 {
                return Err(ReflectionError::InvalidMtz);
            }
            for _ in 0..count {
                history.push(record(bytes, offset)?.trim().into());
                offset += RECORD_BYTES;
            }
        } else if trimmed.starts_with("MTZBATS") {
            return Err(ReflectionError::Unsupported(
                "unmerged MTZ batch headers".into(),
            ));
        }
    }
    Err(ReflectionError::TruncatedMtz)
}
