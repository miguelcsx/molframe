//! A trajectory as one rectangular block of coordinates.

use super::model::TrajectoryData;
use crate::TrajectoryError;

/// Every frame's positions in one contiguous block, frame-major.
#[derive(Clone, Debug, PartialEq)]
pub struct DensePositions {
    /// Number of frames.
    pub frames: usize,
    /// Atoms in every frame.
    pub atoms: usize,
    /// `frames * atoms * 3` coordinates in ångström: frame, then atom, then axis.
    pub coordinates: Vec<f32>,
    /// Simulation time of each frame in picoseconds; `NaN` where the source
    /// records none, never an invented value.
    pub times: Vec<f64>,
}

impl TrajectoryData {
    /// The frames as one rectangular block, for consumers that want a single
    /// array rather than a list of frames.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::AtomCountMismatch`] when frames disagree on
    /// atom count, since they cannot form a rectangle.
    pub fn dense_positions(&self) -> Result<DensePositions, TrajectoryError> {
        let atoms = self.frames.first().map_or(0, |frame| frame.positions.len());
        let mut coordinates = Vec::with_capacity(self.frames.len() * atoms * 3);
        let mut times = Vec::with_capacity(self.frames.len());
        for frame in &self.frames {
            if frame.positions.len() != atoms {
                return Err(TrajectoryError::AtomCountMismatch {
                    expected: atoms,
                    found: frame.positions.len(),
                });
            }
            coordinates.extend(frame.positions.iter().flatten());
            times.push(match frame.time {
                Some(time) => time,
                None => f64::NAN,
            });
        }
        Ok(DensePositions {
            frames: self.frames.len(),
            atoms,
            coordinates,
            times,
        })
    }
}
