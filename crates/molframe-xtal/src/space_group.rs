//! Offline catalogue for all crystallographic Hall settings.

use crate::{Rational, SymmetryOperation, SymmetrySet};
use molframe_core::{Code, Diagnostic};
use std::sync::OnceLock;

const CATALOGUE: &[u8] = include_bytes!("../data/space-groups.bin");
const MAGIC: &[u8; 8] = b"PDBXSG01";
const SETTING_COUNT: usize = 530;

static SETTINGS: OnceLock<Result<Box<[SpaceGroupSetting]>, Diagnostic>> = OnceLock::new();

/// One Hall setting, including its complete exact operator list.
#[derive(Clone, Debug)]
pub struct SpaceGroupSetting {
    /// Catalogue serial number in `1..=530`.
    pub hall_number: u16,
    /// International Tables space-group type in `1..=230`.
    pub international_number: u16,
    /// Short Hermann–Mauguin symbol.
    pub international_short: Box<str>,
    /// Full Hermann–Mauguin symbol.
    pub international_full: Box<str>,
    /// Hall symbol identifying setting and origin uniquely.
    pub hall_symbol: Box<str>,
    /// Unique-axis, origin or cell choice recorded by the catalogue.
    pub choice: Box<str>,
    /// Exact symmetry representatives for this setting.
    pub operations: Box<[SymmetryOperation]>,
}

impl SpaceGroupSetting {
    /// Converts catalogue metadata into structure-attached symmetry data.
    #[must_use]
    pub fn symmetry_set(&self) -> SymmetrySet {
        SymmetrySet {
            hall_number: Some(self.hall_number),
            international_number: Some(self.international_number),
            hermann_mauguin: Some(self.international_short.clone()),
            hall: Some(self.hall_symbol.clone()),
            crystal_system: None,
            choice: (!self.choice.is_empty()).then(|| self.choice.clone()),
            operations: self.operations.to_vec(),
        }
    }
}

/// Looks up one of the 530 settings by its catalogue Hall number.
///
/// # Errors
///
/// Returns `E6018` for a number outside `1..=530` or invalid bundled data.
pub fn space_group_setting(hall_number: u16) -> Result<&'static SpaceGroupSetting, Diagnostic> {
    let index = usize::from(hall_number)
        .checked_sub(1)
        .ok_or_else(catalogue_error)?;
    settings()?.get(index).ok_or_else(catalogue_error)
}

/// Looks up a setting by its unique Hall symbol.
///
/// ASCII whitespace is collapsed before comparison.
///
/// # Errors
///
/// Returns `E6018` when the symbol is absent or bundled data is invalid.
pub fn space_group_by_hall(symbol: &str) -> Result<&'static SpaceGroupSetting, Diagnostic> {
    let requested = normalized(symbol);
    settings()?
        .iter()
        .find(|setting| normalized(&setting.hall_symbol) == requested)
        .ok_or_else(catalogue_error)
}

/// Returns every Hall setting belonging to one International Tables type.
///
/// # Errors
///
/// Returns `E6018` for a type outside `1..=230` or invalid bundled data.
pub fn space_group_settings(
    international_number: u16,
) -> Result<Vec<&'static SpaceGroupSetting>, Diagnostic> {
    if !(1..=230).contains(&international_number) {
        return Err(catalogue_error());
    }
    Ok(settings()?
        .iter()
        .filter(|setting| setting.international_number == international_number)
        .collect())
}

pub(crate) fn resolve_metadata_operations(set: &mut SymmetrySet) -> Result<(), Diagnostic> {
    let setting = if let Some(hall) = set.hall.as_deref() {
        Some(space_group_by_hall(hall)?)
    } else if let Some(number) = set.international_number {
        space_group_settings(number)?.into_iter().next()
    } else {
        None
    };
    let Some(setting) = setting else {
        return Ok(());
    };
    set.hall_number = Some(setting.hall_number);
    set.international_number = Some(setting.international_number);
    if set.hermann_mauguin.is_none() {
        set.hermann_mauguin = Some(setting.international_short.clone());
    }
    if set.hall.is_none() {
        set.hall = Some(setting.hall_symbol.clone());
    }
    set.choice = (!setting.choice.is_empty()).then(|| setting.choice.clone());
    set.operations = setting.operations.to_vec();
    Ok(())
}

fn settings() -> Result<&'static [SpaceGroupSetting], Diagnostic> {
    match SETTINGS.get_or_init(decode_catalogue) {
        Ok(settings) => Ok(settings),
        Err(finding) => Err(finding.clone()),
    }
}

fn decode_catalogue() -> Result<Box<[SpaceGroupSetting]>, Diagnostic> {
    let mut cursor = Cursor::new(CATALOGUE);
    if cursor.take(MAGIC.len())? != MAGIC {
        return Err(catalogue_error());
    }
    if usize::from(cursor.u16()?) != SETTING_COUNT {
        return Err(catalogue_error());
    }
    let mut output = Vec::with_capacity(SETTING_COUNT);
    for expected in 1..=SETTING_COUNT {
        let hall_number = cursor.u16()?;
        if usize::from(hall_number) != expected {
            return Err(catalogue_error());
        }
        let international_number = cursor.u16()?;
        if !(1..=230).contains(&international_number) {
            return Err(catalogue_error());
        }
        let international_short = cursor.text()?;
        let international_full = cursor.text()?;
        let hall_symbol = cursor.text()?;
        let choice = cursor.text()?;
        let operation_count = usize::from(cursor.u16()?);
        let mut operations = Vec::with_capacity(operation_count);
        for operation in 0..operation_count {
            operations.push(cursor.operation(operation)?);
        }
        output.push(SpaceGroupSetting {
            hall_number,
            international_number,
            international_short,
            international_full,
            hall_symbol,
            choice,
            operations: operations.into_boxed_slice(),
        });
    }
    if cursor.remaining() != 0 {
        return Err(catalogue_error());
    }
    Ok(output.into_boxed_slice())
}

struct Cursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Cursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], Diagnostic> {
        let end = self.offset.checked_add(count).ok_or_else(catalogue_error)?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(catalogue_error)?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, Diagnostic> {
        self.take(1)?.first().copied().ok_or_else(catalogue_error)
    }

    fn i8(&mut self) -> Result<i8, Diagnostic> {
        Ok(i8::from_le_bytes([self.u8()?]))
    }

    fn u16(&mut self) -> Result<u16, Diagnostic> {
        let bytes: [u8; 2] = self.take(2)?.try_into().map_err(|_| catalogue_error())?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn text(&mut self) -> Result<Box<str>, Diagnostic> {
        let length = usize::from(self.u8()?);
        let value = std::str::from_utf8(self.take(length)?).map_err(|_| catalogue_error())?;
        Ok(value.into())
    }

    fn operation(&mut self, index: usize) -> Result<SymmetryOperation, Diagnostic> {
        let mut rotation = [[0; 3]; 3];
        for row in &mut rotation {
            for value in row {
                *value = i32::from(self.i8()?);
            }
        }
        let mut translation = [Rational::ZERO; 3];
        for value in &mut translation {
            *value = Rational::new(i32::from(self.u8()?), 12)?;
        }
        Ok(SymmetryOperation {
            id: (index + 1).to_string().into_boxed_str(),
            rotation,
            translation,
        })
    }

    fn remaining(&self) -> usize {
        if self.offset >= self.bytes.len() {
            0
        } else {
            self.bytes.len() - self.offset
        }
    }
}

fn normalized(value: &str) -> String {
    value.split_ascii_whitespace().collect::<Vec<_>>().join(" ")
}

fn catalogue_error() -> Diagnostic {
    Diagnostic::new(Code::E6018)
}

#[cfg(test)]
#[path = "space_group_tests.rs"]
mod tests;
