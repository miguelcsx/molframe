//! Shannon ionic radii indexed by element, charge and coordination.

use molframe_core::{Code, Diagnostic, Element};
use std::ops::Range;
use std::sync::OnceLock;

const DATA: &[u8] = include_bytes!("../data/ionic-radii.bin");
const MAGIC: &[u8; 8] = b"PDBXIR01";
static CATALOGUE: OnceLock<Result<IonicCatalogue, Diagnostic>> = OnceLock::new();

/// Spin state attached to a Shannon ionic-radius entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IonicSpin {
    /// No spin distinction applies or was supplied.
    Unspecified,
    /// High-spin configuration.
    High,
    /// Low-spin configuration.
    Low,
}

/// One charge- and coordination-specific ionic radius.
#[derive(Clone, Debug, PartialEq)]
pub struct IonicRadius {
    /// Formal ionic charge.
    pub charge: i8,
    /// Coordination geometry/count in Shannon's notation (`VI`, `IVSQ`, …).
    pub coordination: Box<str>,
    /// Spin state where the table distinguishes it.
    pub spin: IonicSpin,
    /// Effective ionic radius in ångström.
    pub ionic_radius: f32,
    /// Crystal radius in ångström.
    pub crystal_radius: f32,
    /// Source marks this entry as among the most reliable, when specified.
    pub most_reliable: Option<bool>,
}

/// Returns all Shannon entries for one element in deterministic table order.
///
/// Missing entries produce an empty slice rather than a substituted estimate.
///
/// # Errors
///
/// Returns `E3201` when the bundled reference table is invalid.
pub fn ionic_radii(element: Element) -> Result<&'static [IonicRadius], Diagnostic> {
    let catalogue = catalogue()?;
    let range = catalogue
        .ranges
        .get(element.atomic_number() as usize)
        .cloned()
        .ok_or_else(data_error)?;
    catalogue.records.get(range).ok_or_else(data_error)
}

struct IonicCatalogue {
    records: Box<[IonicRadius]>,
    ranges: Box<[Range<usize>]>,
}

fn catalogue() -> Result<&'static IonicCatalogue, Diagnostic> {
    match CATALOGUE.get_or_init(decode) {
        Ok(catalogue) => Ok(catalogue),
        Err(finding) => Err(finding.clone()),
    }
}

fn decode() -> Result<IonicCatalogue, Diagnostic> {
    let mut cursor = Cursor::new(DATA);
    if cursor.take(MAGIC.len())? != MAGIC {
        return Err(data_error());
    }
    let count = usize::from(cursor.u16()?);
    let mut records = Vec::with_capacity(count);
    let mut starts = [usize::MAX; 119];
    let mut ends = [0usize; 119];
    let mut previous = 0u8;
    for _ in 0..count {
        let atomic_number = cursor.u8()?;
        if atomic_number == 0 || atomic_number > 118 || atomic_number < previous {
            return Err(data_error());
        }
        previous = atomic_number;
        let index = usize::from(atomic_number);
        starts[index] = starts[index].min(records.len());
        records.push(IonicRadius {
            charge: cursor.i8()?,
            coordination: cursor.text()?,
            spin: spin(&cursor.text()?)?,
            ionic_radius: cursor.f32()?,
            crystal_radius: cursor.f32()?,
            most_reliable: reliability(cursor.u8()?)?,
        });
        ends[index] = records.len();
    }
    if cursor.remaining() != 0 {
        return Err(data_error());
    }
    let ranges = starts
        .into_iter()
        .zip(ends)
        .map(|(start, end)| {
            if start == usize::MAX {
                0..0
            } else {
                start..end
            }
        })
        .collect();
    Ok(IonicCatalogue {
        records: records.into_boxed_slice(),
        ranges,
    })
}

fn spin(value: &str) -> Result<IonicSpin, Diagnostic> {
    match value {
        "" => Ok(IonicSpin::Unspecified),
        "HS" => Ok(IonicSpin::High),
        "LS" => Ok(IonicSpin::Low),
        _ => Err(data_error()),
    }
}

fn reliability(value: u8) -> Result<Option<bool>, Diagnostic> {
    match value {
        0 => Ok(Some(false)),
        1 => Ok(Some(true)),
        2 => Ok(None),
        _ => Err(data_error()),
    }
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
        let end = self.offset.checked_add(count).ok_or_else(data_error)?;
        let value = self.bytes.get(self.offset..end).ok_or_else(data_error)?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, Diagnostic> {
        self.take(1)?.first().copied().ok_or_else(data_error)
    }

    fn i8(&mut self) -> Result<i8, Diagnostic> {
        Ok(i8::from_le_bytes([self.u8()?]))
    }

    fn u16(&mut self) -> Result<u16, Diagnostic> {
        let value: [u8; 2] = self.take(2)?.try_into().map_err(|_| data_error())?;
        Ok(u16::from_le_bytes(value))
    }

    fn f32(&mut self) -> Result<f32, Diagnostic> {
        let value: [u8; 4] = self.take(4)?.try_into().map_err(|_| data_error())?;
        let value = f32::from_le_bytes(value);
        value.is_finite().then_some(value).ok_or_else(data_error)
    }

    fn text(&mut self) -> Result<Box<str>, Diagnostic> {
        let length = usize::from(self.u8()?);
        let value = std::str::from_utf8(self.take(length)?).map_err(|_| data_error())?;
        Ok(value.into())
    }

    fn remaining(&self) -> usize {
        if self.offset >= self.bytes.len() {
            0
        } else {
            self.bytes.len() - self.offset
        }
    }
}

fn data_error() -> Diagnostic {
    Diagnostic::new(Code::E3201).with_context("reference", "ionic radii")
}

#[cfg(test)]
#[path = "ionic_tests.rs"]
mod tests;
