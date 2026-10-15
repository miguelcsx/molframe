//! Offline catalogue for all crystallographic Hall settings.
//!
//! The 530 settings are definitional: the catalogue of space-group types and
//! their exact operators is fixed by crystallography, so no profile can swap
//! it, and decoding it from a missing file would make every symmetry query
//! depend on a download. It therefore ships inside the binary. Cost: one
//! embedded blob decoded once on first use, then shared for the process.

use crate::{Rational, SymmetryOperation, SymmetrySet};
use molframe_core::{Code, Diagnostic};
use std::sync::OnceLock;

#[path = "space_group_resolve.rs"]
mod resolve;

pub(crate) use resolve::resolve_metadata_operations;
pub use resolve::space_group_by_hermann_mauguin_in_cell;

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

/// Looks up a setting by its Hall symbol.
///
/// ASCII whitespace is collapsed before comparison. A change-of-basis suffix,
/// as an algebraic triplet (`P 61 2 (x,y,z+5/12)`) or as shifts in twelfths
/// (`P 61 2 (0 0 5)`), is applied to the base symbol and the result found in
/// the catalogue whichever way the file spells it.
///
/// # Errors
///
/// Returns `E6018` when the symbol is absent, its change of basis leads to a
/// setting outside the catalogue, or bundled data is invalid.
pub fn space_group_by_hall(symbol: &str) -> Result<&'static SpaceGroupSetting, Diagnostic> {
    match resolve::resolve_hall(symbol)? {
        resolve::HallResolution::Setting(setting) => Ok(setting),
        resolve::HallResolution::Derived { .. } => Err(catalogue_error()),
    }
}

/// Looks up a setting by its Hermann–Mauguin symbol, short or full.
///
/// Case, whitespace and screw-axis underscores are ignored, so the spelling of
/// a PDB `CRYST1` record (`P 21 21 21`) resolves like the catalogue's. A
/// `:choice` suffix (`R 3 :H`, `P 4/n :2`) selects among settings sharing the
/// symbol, and the PDB's `H` lattice letter (`H 3 2`) is the rhombohedral type
/// on hexagonal axes. Otherwise, when a symbol belongs to several settings the
/// first catalogue entry, the standard setting of that type, is returned; use
/// [`space_group_by_hermann_mauguin_in_cell`] to let the cell decide the axes
/// of a rhombohedral group.
///
/// # Errors
///
/// Returns `E6018` when no setting carries the symbol or bundled data is invalid.
pub fn space_group_by_hermann_mauguin(
    symbol: &str,
) -> Result<&'static SpaceGroupSetting, Diagnostic> {
    let key = symbol_key(symbol);
    let (name, mut choice) = match key.split_once(':') {
        Some((name, choice)) => (name.to_owned(), Some(choice.to_owned())),
        None => (key, None),
    };
    let name = match name.strip_prefix('H') {
        Some(rest) => {
            choice.get_or_insert_with(|| "H".to_owned());
            format!("R{rest}")
        }
        None => name,
    };
    settings()?
        .iter()
        .find(|setting| {
            (symbol_key(&setting.international_full) == name
                || symbol_key(&setting.international_short) == name)
                && choice
                    .as_deref()
                    .is_none_or(|wanted| setting.choice.to_ascii_uppercase() == wanted)
        })
        .ok_or_else(catalogue_error)
}

/// A symbol with case, spacing and screw-axis subscripts removed, so the
/// catalogue's `P 2_1 2_1 2_1` and a PDB file's `P 21 21 21` compare equal.
fn symbol_key(symbol: &str) -> String {
    symbol
        .chars()
        .filter(|character| !character.is_whitespace() && *character != '_')
        .map(|character| character.to_ascii_uppercase())
        .collect()
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
