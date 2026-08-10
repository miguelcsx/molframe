//! GROMACS TRR trajectory reader using the format's native XDR encoding.

use crate::Timestep;
use crate::cell::cell_from_vectors;
use crate::numeric::f32_from_f64;

const MAGIC: i32 = 1993;
const VERSION: &str = "GMX_trn_file";
const NM_TO_ANGSTROM: f32 = 10.0;
const FORCE_TO_CANONICAL: f32 = 0.1;

/// Floating-point representation carried by a TRR frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrrPrecision {
    /// IEEE-754 binary32 values.
    Single,
    /// IEEE-754 binary64 values.
    Double,
}

impl TrrPrecision {
    pub(crate) const fn bytes(self) -> usize {
        match self {
            Self::Single => 4,
            Self::Double => 8,
        }
    }
}

/// Parsed TRR trajectory and the precision of every frame.
#[derive(Clone, Debug, PartialEq)]
pub struct TrrTrajectory {
    /// Frames in file order.
    pub frames: Vec<Timestep>,
    /// Simulation step stored in each frame.
    pub steps: Vec<i32>,
    /// Per-frame precision, since concatenated TRR streams may vary it.
    pub precision: Vec<TrrPrecision>,
}

/// Malformed or unsupported TRR data.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum TrrError {
    /// Input ended inside an XDR value or declared block.
    #[error("truncated TRR data at byte {offset}")]
    Truncated {
        /// First unavailable byte.
        offset: usize,
    },
    /// Magic number or version string is not a supported TRR header.
    #[error("invalid TRR header")]
    InvalidHeader,
    /// Counts, precision or block sizes contradict one another.
    #[error("inconsistent TRR block sizes")]
    InvalidSizes,
    /// A frame contains a non-finite value or invalid unit cell.
    #[error("invalid numeric value in TRR frame")]
    InvalidValue,
}

/// Parses all concatenated frames from a TRR stream.
///
/// Positions and cells are converted from nm to ångström, velocities from
/// nm/ps to ångström/ps, and forces from kJ mol⁻¹ nm⁻¹ to kJ mol⁻¹ Å⁻¹.
///
/// # Errors
///
/// Returns the first header, XDR, size or numeric error without publishing a
/// partial trajectory.
pub fn parse_trr(bytes: &[u8]) -> Result<TrrTrajectory, TrrError> {
    let mut reader = Xdr::new(bytes);
    let mut frames = Vec::new();
    let mut steps = Vec::new();
    let mut precisions = Vec::new();
    while !reader.is_empty() {
        let header = read_header(&mut reader)?;
        let precision = header.precision()?;
        skip_prefix_blocks(&mut reader, &header)?;
        let cell = read_cell(&mut reader, &header, precision)?;
        reader.skip(header.virial() + header.pressure() + header.topology() + header.symmetry())?;
        let positions = read_vectors(
            &mut reader,
            header.positions(),
            header.atoms,
            precision,
            NM_TO_ANGSTROM,
        )?
        .ok_or(TrrError::InvalidSizes)?;
        let velocities = read_vectors(
            &mut reader,
            header.velocities(),
            header.atoms,
            precision,
            NM_TO_ANGSTROM,
        )?;
        let forces = read_vectors(
            &mut reader,
            header.forces(),
            header.atoms,
            precision,
            FORCE_TO_CANONICAL,
        )?;
        let frame = frames.len();
        frames.push(Timestep {
            frame,
            time: Some(header.time),
            positions,
            velocities,
            forces,
            cell,
            ..Timestep::default()
        });
        steps.push(header.step);
        precisions.push(precision);
    }
    populate_dt(&mut frames);
    Ok(TrrTrajectory {
        frames,
        steps,
        precision: precisions,
    })
}

struct Header {
    sizes: [usize; 10],
    atoms: usize,
    step: i32,
    time: f64,
}

impl Header {
    fn precision(&self) -> Result<TrrPrecision, TrrError> {
        let scalar_bytes = if self.sizes[2] != 0 {
            exact_div(self.sizes[2], 9)?
        } else {
            let vector_size = self.sizes[7..]
                .iter()
                .copied()
                .find(|size| *size != 0)
                .ok_or(TrrError::InvalidSizes)?;
            exact_div(
                vector_size,
                self.atoms.checked_mul(3).ok_or(TrrError::InvalidSizes)?,
            )?
        };
        match scalar_bytes {
            4 => Ok(TrrPrecision::Single),
            8 => Ok(TrrPrecision::Double),
            _ => Err(TrrError::InvalidSizes),
        }
    }

    const fn input_record(&self) -> usize {
        self.sizes[0]
    }
    const fn energies(&self) -> usize {
        self.sizes[1]
    }
    const fn cell(&self) -> usize {
        self.sizes[2]
    }
    const fn virial(&self) -> usize {
        self.sizes[3]
    }
    const fn pressure(&self) -> usize {
        self.sizes[4]
    }
    const fn topology(&self) -> usize {
        self.sizes[5]
    }
    const fn symmetry(&self) -> usize {
        self.sizes[6]
    }
    const fn positions(&self) -> usize {
        self.sizes[7]
    }
    const fn velocities(&self) -> usize {
        self.sizes[8]
    }
    const fn forces(&self) -> usize {
        self.sizes[9]
    }
}

fn read_header(reader: &mut Xdr<'_>) -> Result<Header, TrrError> {
    if reader.i32()? != MAGIC {
        return Err(TrrError::InvalidHeader);
    }
    let declared_version_length =
        usize::try_from(reader.i32()?).map_err(|_| TrrError::InvalidHeader)?;
    if declared_version_length != VERSION.len() + 1 || reader.string()? != VERSION {
        return Err(TrrError::InvalidHeader);
    }
    let mut sizes = [0usize; 10];
    for size in &mut sizes {
        *size = usize::try_from(reader.i32()?).map_err(|_| TrrError::InvalidSizes)?;
    }
    let atoms = usize::try_from(reader.i32()?).map_err(|_| TrrError::InvalidSizes)?;
    let step = reader.i32()?;
    let _integration_records = reader.i32()?;
    if atoms == 0 || sizes[7] == 0 {
        return Err(TrrError::InvalidSizes);
    }
    let precision = Header {
        sizes,
        atoms,
        step,
        time: 0.0,
    }
    .precision()?;
    let time = reader.real(precision)?;
    let lambda = reader.real(precision)?;
    if !time.is_finite() || !lambda.is_finite() {
        return Err(TrrError::InvalidValue);
    }
    Ok(Header {
        sizes,
        atoms,
        step,
        time,
    })
}

fn skip_prefix_blocks(reader: &mut Xdr<'_>, header: &Header) -> Result<(), TrrError> {
    reader.skip(header.input_record() + header.energies())
}

fn read_cell(
    reader: &mut Xdr<'_>,
    header: &Header,
    precision: TrrPrecision,
) -> Result<Option<pdbiox_core::structure::UnitCell>, TrrError> {
    if header.cell() == 0 {
        return Ok(None);
    }
    if header.cell() != 9 * precision.bytes() {
        return Err(TrrError::InvalidSizes);
    }
    let mut vectors = [[0.0; 3]; 3];
    for row in &mut vectors {
        for value in row {
            *value = reader.real(precision)? * f64::from(NM_TO_ANGSTROM);
        }
    }
    cell_from_vectors(vectors)
        .ok_or(TrrError::InvalidValue)
        .map(Some)
}

fn read_vectors(
    reader: &mut Xdr<'_>,
    size: usize,
    atoms: usize,
    precision: TrrPrecision,
    scale: f32,
) -> Result<Option<Vec<[f32; 3]>>, TrrError> {
    if size == 0 {
        return Ok(None);
    }
    if size != atoms * 3 * precision.bytes() {
        return Err(TrrError::InvalidSizes);
    }
    let mut vectors = Vec::with_capacity(atoms);
    for _ in 0..atoms {
        let vector = [
            f32_from_f64(reader.real(precision)?).ok_or(TrrError::InvalidValue)? * scale,
            f32_from_f64(reader.real(precision)?).ok_or(TrrError::InvalidValue)? * scale,
            f32_from_f64(reader.real(precision)?).ok_or(TrrError::InvalidValue)? * scale,
        ];
        if vector.iter().any(|value| !value.is_finite()) {
            return Err(TrrError::InvalidValue);
        }
        vectors.push(vector);
    }
    Ok(Some(vectors))
}

fn populate_dt(frames: &mut [Timestep]) {
    for index in 1..frames.len() {
        if let (Some(previous), Some(current)) = (frames[index - 1].time, frames[index].time) {
            frames[index].dt = Some(current - previous);
        }
    }
    if frames.len() > 1 {
        frames[0].dt = frames[1].dt;
    }
}

fn exact_div(value: usize, divisor: usize) -> Result<usize, TrrError> {
    if divisor == 0 || !value.is_multiple_of(divisor) {
        Err(TrrError::InvalidSizes)
    } else {
        Ok(value / divisor)
    }
}

struct Xdr<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Xdr<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }
    fn is_empty(&self) -> bool {
        self.cursor == self.bytes.len()
    }
    fn take(&mut self, length: usize) -> Result<&'a [u8], TrrError> {
        let start = self.cursor;
        let end = start
            .checked_add(length)
            .ok_or(TrrError::Truncated { offset: start })?;
        let value = self
            .bytes
            .get(start..end)
            .ok_or(TrrError::Truncated { offset: start })?;
        self.cursor = end;
        Ok(value)
    }
    fn skip(&mut self, length: usize) -> Result<(), TrrError> {
        self.take(length).map(|_| ())
    }
    fn i32(&mut self) -> Result<i32, TrrError> {
        let offset = self.cursor;
        let raw = self
            .take(4)?
            .try_into()
            .map_err(|_| TrrError::Truncated { offset })?;
        Ok(i32::from_be_bytes(raw))
    }
    fn real(&mut self, precision: TrrPrecision) -> Result<f64, TrrError> {
        let offset = self.cursor;
        match precision {
            TrrPrecision::Single => {
                let raw = self
                    .take(4)?
                    .try_into()
                    .map_err(|_| TrrError::Truncated { offset })?;
                Ok(f64::from(f32::from_be_bytes(raw)))
            }
            TrrPrecision::Double => {
                let raw = self
                    .take(8)?
                    .try_into()
                    .map_err(|_| TrrError::Truncated { offset })?;
                Ok(f64::from_be_bytes(raw))
            }
        }
    }
    fn string(&mut self) -> Result<&'a str, TrrError> {
        let length = usize::try_from(self.i32()?).map_err(|_| TrrError::InvalidHeader)?;
        let bytes = self.take(length)?;
        let padding = (4 - length % 4) % 4;
        self.skip(padding)?;
        std::str::from_utf8(bytes).map_err(|_| TrrError::InvalidHeader)
    }
}

#[cfg(test)]
#[path = "trr_tests.rs"]
mod tests;
