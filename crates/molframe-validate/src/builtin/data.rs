//! Decoders and lookups for Ramachandran and rotamer density tables.
//!
//! Generic over the table contents: the decoders accept the binary formats
//! documented in `data/README.md` from any source, bundled or user supplied,
//! and report damage as an error value rather than trusting their input. The
//! bundled bytes live in the `bundled` module, not here.

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

    /// Density below which a residue is an outlier.
    #[must_use]
    pub const fn outlier_density(self) -> f32 {
        match self {
            Self::General => 0.0005,
            Self::CisPro => 0.002,
            _ => 0.001,
        }
    }
}

/// Density at or above which a Ramachandran residue is favoured.
pub const RAMA_FAVOURED: f32 = 0.02;

/// Density cut-offs that turn a Ramachandran density into a class.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RamaThresholds {
    /// Density at or above which a residue is favoured.
    pub favoured: f32,
    /// Density below which a residue is an outlier, per category code.
    pub outlier: [f32; 6],
}

impl Default for RamaThresholds {
    fn default() -> Self {
        let outlier = [
            RamaCategory::General,
            RamaCategory::Gly,
            RamaCategory::TransPro,
            RamaCategory::CisPro,
            RamaCategory::PrePro,
            RamaCategory::IleVal,
        ]
        .map(RamaCategory::outlier_density);
        Self {
            favoured: RAMA_FAVOURED,
            outlier,
        }
    }
}

/// Six decoded Ramachandran grids, as densities, with their cut-offs.
#[derive(Debug)]
pub struct RamaGrids {
    grids: Vec<Vec<f32>>,
    thresholds: RamaThresholds,
}

impl RamaGrids {
    /// Decodes a Ramachandran table file.
    ///
    /// # Errors
    ///
    /// Returns [`TableError`] if the bytes are not a well-formed table.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, TableError> {
        decode_rama(bytes)
    }

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

/// Failure to decode a reference table.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[error("reference table `{0}` is malformed")]
pub struct TableError(pub &'static str);

fn dequantise(code: u8) -> f32 {
    if code == 0 {
        return 0.0;
    }
    let exponent = 6.0 * (f64::from(code) - 1.0) / 254.0 - 6.0;
    crate::numeric::f64_to_f32(10f64.powf(exponent))
}

fn decode_rama(bytes: &[u8]) -> Result<RamaGrids, TableError> {
    let bad = TableError("rama");
    if bytes.get(..4) != Some(b"MFRM") {
        return Err(bad);
    }
    let count = usize::from(*bytes.get(4).ok_or(bad)?);
    let cells = RAMA_BINS * RAMA_BINS;
    let mut grids = Vec::with_capacity(count);
    let mut at = 5;
    for _ in 0..count {
        let body = bytes.get(at + 8..at + 8 + cells).ok_or(bad)?;
        grids.push(body.iter().map(|code| dequantise(*code)).collect());
        at += 8 + cells;
    }
    if grids.len() < 6 {
        return Err(bad);
    }
    Ok(RamaGrids {
        grids,
        thresholds: RamaThresholds::default(),
    })
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
    tables: Vec<(String, RotamerTable)>,
}

impl RotamerTables {
    /// Decodes a rotamer table file.
    ///
    /// # Errors
    ///
    /// Returns [`TableError`] if the bytes are not a well-formed table.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, TableError> {
        decode_rotamers(bytes)
    }

    /// The table with this lower-case name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&RotamerTable> {
        self.tables
            .iter()
            .find(|(label, _)| label == name)
            .map(|(_, table)| table)
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

fn decode_rotamers(bytes: &[u8]) -> Result<RotamerTables, TableError> {
    let bad = TableError("rotamer");
    if bytes.get(..4) != Some(b"MFRT") {
        return Err(bad);
    }
    let mut cursor = Cursor { bytes, at: 4 };
    let count = cursor.byte().ok_or(bad)?;
    let mut tables = Vec::with_capacity(usize::from(count));
    for _ in 0..count {
        tables.push(decode_table(&mut cursor).ok_or(bad)?);
    }
    Ok(RotamerTables { tables })
}

#[cfg(test)]
#[path = "data_tests.rs"]
mod tests;
