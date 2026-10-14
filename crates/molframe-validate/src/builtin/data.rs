//! Decoders and lookups for Ramachandran and rotamer density tables.
//!
//! Generic over the table contents: the decoders accept the binary formats
//! documented beside the data from any source and report damage as an error
//! value rather than trusting their input. Each file starts with the shared
//! identity header (name, version, source, licence), so a result can say which
//! table produced it. A Ramachandran header also carries its class cut-offs
//! (`favoured`, and `outlier` per category code). The library compiles in no
//! table and no cut-off.
//!
//! Complexity: decoding is O(table bytes) once per table, then shared through an
//! `Arc`; a density or rotamer lookup is O(dimensions) with direct grid indexing
//! and no allocation. A rotamer table is found by name in O(1) through a hash
//! index.

use std::collections::HashMap;

use molframe_chem::{
    IdentityOnly, Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, split_binary,
};

const RAMA_MAGIC: [u8; 4] = *b"MFR2";
const ROTAMER_MAGIC: [u8; 4] = *b"MFT2";

/// Cells per axis of a Ramachandran grid.
pub const RAMA_BINS: usize = 180;
const RAMA_CELL: f64 = 2.0;

/// Which Ramachandran distribution a residue is scored against.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RamaCategory {
    /// Every residue not covered below.
    General,
    /// Glycine.
    Gly,
    /// Proline preceded by a trans peptide.
    TransPro,
    /// Proline preceded by a cis peptide.
    CisPro,
    /// Any residue followed by proline, except Gly and Pro.
    PrePro,
    /// Isoleucine and valine.
    IleVal,
}

impl RamaCategory {
    /// Stable code used in the per-residue arrays.
    #[must_use]
    pub const fn code(self) -> u8 {
        match self {
            Self::General => 0,
            Self::Gly => 1,
            Self::TransPro => 2,
            Self::CisPro => 3,
            Self::PrePro => 4,
            Self::IleVal => 5,
        }
    }
}

/// Density cut-offs that turn a Ramachandran density into a class.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RamaThresholds {
    /// Density at or above which a residue is favoured.
    pub favoured: f32,
    /// Density below which a residue is an outlier, per category code.
    pub outlier: [f32; 6],
}

/// Reads the `favoured` and `outlier` headers of a Ramachandran table.
#[derive(Default)]
struct ThresholdHeader {
    favoured: Option<f32>,
    outlier: Option<[f32; 6]>,
}

const RAMA: &str = "ramachandran table";

fn cutoff(line: usize, field: &str) -> Result<f32, TableTextError> {
    let row = Row::new(RAMA, line);
    let value = row.finite::<f32>(field)?;
    if value < 0.0 {
        return Err(row.fail(format!("`{field}` is not a non-negative number")));
    }
    Ok(value)
}

impl TableLines for ThresholdHeader {
    fn header(&mut self, line: usize, key: &str, value: &str) -> Result<bool, TableTextError> {
        match key {
            "favoured" => self.favoured = Some(cutoff(line, value.trim())?),
            "outlier" => {
                let mut cuts = [0.0f32; 6];
                let mut fields = value.split_whitespace();
                for cut in &mut cuts {
                    let Some(field) = fields.next() else {
                        return Err(TableTextError::at(
                            "ramachandran table",
                            line,
                            "`outlier` needs six cut-offs, one per category code",
                        ));
                    };
                    *cut = cutoff(line, field)?;
                }
                if fields.next().is_some() {
                    return Err(TableTextError::at(
                        "ramachandran table",
                        line,
                        "`outlier` has more than six cut-offs",
                    ));
                }
                self.outlier = Some(cuts);
            }
            _ => return Ok(false),
        }
        Ok(true)
    }

    fn row(&mut self, line: usize, _: &str) -> Result<(), TableTextError> {
        Err(TableTextError::at(
            "ramachandran table",
            line,
            "the header takes no rows",
        ))
    }
}

/// Six decoded Ramachandran grids, as densities, with their cut-offs.
#[derive(Debug)]
pub struct RamaGrids {
    grids: Vec<Vec<f32>>,
    thresholds: RamaThresholds,
}

impl RamaGrids {
    /// Replaces the class cut-offs.
    #[must_use]
    pub fn with_thresholds(mut self, thresholds: RamaThresholds) -> Self {
        self.thresholds = thresholds;
        self
    }

    /// The class cut-offs in force.
    #[must_use]
    pub const fn thresholds(&self) -> RamaThresholds {
        self.thresholds
    }

    /// Class code (0 outlier, 1 allowed, 2 favoured) for a density.
    #[must_use]
    pub fn classify(&self, category: RamaCategory, density: f32) -> u8 {
        let outlier = match self.thresholds.outlier.get(usize::from(category.code())) {
            Some(cut) => *cut,
            None => 0.0,
        };
        if density < outlier {
            0
        } else if density >= self.thresholds.favoured {
            2
        } else {
            1
        }
    }

    /// The density at (phi, psi) in degrees for one category.
    #[must_use]
    pub fn density(&self, category: RamaCategory, phi: f64, psi: f64) -> f32 {
        let Some(grid) = self.grids.get(usize::from(category.code())) else {
            return 0.0;
        };
        let cell = bin(phi) * RAMA_BINS + bin(psi);
        match grid.get(cell) {
            Some(density) => *density,
            None => 0.0,
        }
    }
}

fn bin(angle: f64) -> usize {
    let wrapped = (angle + 180.0).rem_euclid(360.0);
    let cell = (wrapped / RAMA_CELL).floor();
    if cell.is_finite() && cell >= 0.0 {
        crate::numeric::f64_to_usize(cell).min(RAMA_BINS - 1)
    } else {
        0
    }
}

fn dequantise(code: u8) -> f32 {
    if code == 0 {
        return 0.0;
    }
    let exponent = 6.0 * (f64::from(code) - 1.0) / 254.0 - 6.0;
    crate::numeric::f64_to_f32(10f64.powf(exponent))
}

impl TableKind for RamaGrids {
    const KIND: &'static str = "ramachandran";
    const EXTENSION: &'static str = "rama";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let what = "ramachandran table";
        let mut header = ThresholdHeader::default();
        let (identity, payload) = split_binary(what, RAMA_MAGIC, bytes, &mut header)?;
        let bad = |message: &str| TableTextError::at(what, 0, message);
        let (Some(favoured), Some(outlier)) = (header.favoured, header.outlier) else {
            return Err(bad(
                "the header must give `favoured` and `outlier` cut-offs",
            ));
        };
        let count = usize::from(*payload.first().ok_or_else(|| bad("no grid count"))?);
        let cells = RAMA_BINS * RAMA_BINS;
        let mut grids = Vec::with_capacity(count);
        let mut at = 1;
        for _ in 0..count {
            let body = payload
                .get(at + 8..at + 8 + cells)
                .ok_or_else(|| bad("a grid is truncated"))?;
            grids.push(body.iter().map(|code| dequantise(*code)).collect());
            at += 8 + cells;
        }
        if grids.len() < 6 {
            return Err(bad("fewer than six grids"));
        }
        let table = Self {
            grids,
            thresholds: RamaThresholds { favoured, outlier },
        };
        Ok(Provided { identity, table })
    }
}

/// Classification of a side-chain conformation cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotamerClass {
    /// Below the allowed contour.
    Outlier,
    /// Between the allowed and favoured contours.
    Allowed,
    /// At or above the favoured contour.
    Favoured,
}

impl RotamerClass {
    /// Stable code used in the per-residue arrays.
    #[must_use]
    pub const fn code(self) -> u8 {
        match self {
            Self::Outlier => 0,
            Self::Allowed => 1,
            Self::Favoured => 2,
        }
    }
}

/// One chi-space classification table.
#[derive(Debug)]
pub struct RotamerTable {
    dims: Vec<(usize, f64)>,
    cells: Vec<u8>,
}

impl RotamerTable {
    /// Number of chi angles this table is indexed by.
    #[must_use]
    pub fn dimensions(&self) -> usize {
        self.dims.len()
    }

    /// Classifies a conformation given by chi angles in degrees.
    ///
    /// Returns `None` when fewer angles are given than the table needs.
    #[must_use]
    pub fn classify(&self, chi: &[f64]) -> Option<RotamerClass> {
        if chi.len() < self.dims.len() {
            return None;
        }
        let mut index = 0usize;
        for ((bins, span), angle) in self.dims.iter().zip(chi) {
            let wrapped = angle.rem_euclid(*span);
            let width = span / crate::numeric::usize_to_f64(*bins);
            let cell = crate::numeric::f64_to_usize((wrapped / width).floor()).min(bins - 1);
            index = index * bins + cell;
        }
        Some(match self.cells.get(index)? {
            2 => RotamerClass::Favoured,
            1 => RotamerClass::Allowed,
            _ => RotamerClass::Outlier,
        })
    }
}

/// Every rotamer table, looked up by name.
#[derive(Debug)]
pub struct RotamerTables {
    tables: HashMap<String, RotamerTable>,
}

impl RotamerTables {
    /// The table with this lower-case name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&RotamerTable> {
        self.tables.get(name)
    }

    /// Number of tables.
    #[must_use]
    pub fn len(&self) -> usize {
        self.tables.len()
    }

    /// Whether there are no tables.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.tables.is_empty()
    }
}

struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, count: usize) -> Option<&'a [u8]> {
        let slice = self.bytes.get(self.at..self.at.checked_add(count)?)?;
        self.at += count;
        Some(slice)
    }

    fn byte(&mut self) -> Option<u8> {
        self.take(1).and_then(|slice| slice.first().copied())
    }

    fn varint(&mut self) -> Option<usize> {
        let mut value = 0usize;
        let mut shift = 0u32;
        loop {
            let byte = self.byte()?;
            value |= usize::from(byte & 0x7F).checked_shl(shift)?;
            if byte & 0x80 == 0 {
                return Some(value);
            }
            shift += 7;
        }
    }
}

fn decode_table(cursor: &mut Cursor<'_>) -> Option<(String, RotamerTable)> {
    let length = usize::from(cursor.byte()?);
    let name = String::from_utf8(cursor.take(length)?.to_vec()).ok()?;
    let dimensions = usize::from(cursor.byte()?);
    let mut dims = Vec::with_capacity(dimensions);
    let mut total = 1usize;
    for _ in 0..dimensions {
        let head = cursor.take(3)?;
        let bins = usize::from(u16::from_le_bytes([*head.first()?, *head.get(1)?]));
        let half = *head.get(2)? == 1;
        total = total.checked_mul(bins)?;
        dims.push((bins, if half { 180.0 } else { 360.0 }));
    }
    let size = cursor.take(4)?;
    let size = u32::from_le_bytes([*size.first()?, *size.get(1)?, *size.get(2)?, *size.get(3)?]);
    let payload = cursor.take(usize::try_from(size).ok()?)?;
    let mut inner = Cursor {
        bytes: payload,
        at: 0,
    };
    let mut cells = Vec::with_capacity(total);
    while inner.at < payload.len() {
        let word = inner.varint()?;
        let class = u8::try_from(word & 3).ok()?;
        cells.resize(cells.len().checked_add(word >> 2)?, class);
    }
    (cells.len() == total).then_some((name, RotamerTable { dims, cells }))
}

impl TableKind for RotamerTables {
    const KIND: &'static str = "rotamer";
    const EXTENSION: &'static str = "rotamer";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let what = "rotamer table";
        let (identity, payload) = split_binary(what, ROTAMER_MAGIC, bytes, &mut IdentityOnly)?;
        let bad = |message: &str| TableTextError::at(what, 0, message);
        let mut cursor = Cursor {
            bytes: payload,
            at: 0,
        };
        let count = cursor.byte().ok_or_else(|| bad("no table count"))?;
        let mut tables = HashMap::with_capacity(usize::from(count));
        for _ in 0..count {
            let (name, table) =
                decode_table(&mut cursor).ok_or_else(|| bad("a table is damaged"))?;
            if tables.insert(name, table).is_some() {
                return Err(bad("a table name is repeated"));
            }
        }
        Ok(Provided {
            identity,
            table: Self { tables },
        })
    }
}

#[cfg(test)]
#[path = "data_tests.rs"]
mod tests;
