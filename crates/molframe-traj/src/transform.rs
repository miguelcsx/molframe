//! Pure on-the-fly coordinate transforms composed over a reader.

use crate::reader::{RandomAccess, Timestep, TrajectoryError, TrajectoryReader, Units};
use molframe_geom::Rigid;
use std::fmt;

/// One pure in-place operation applied after a frame is read.
pub trait FrameTransform: fmt::Debug + Send + Sync {
    /// Stable name recorded by the pipeline.
    fn name(&self) -> &'static str;
    /// Applies the operation to one timestep.
    ///
    /// # Errors
    ///
    /// Returns a selection or degenerate-geometry error without changing
    /// reader position.
    fn apply(&self, timestep: &mut Timestep) -> Result<(), TrajectoryError>;
}

/// Translation or rotation represented by a shared geometric kernel.
#[derive(Clone, Copy, Debug)]
pub struct RigidTransform {
    transform: Rigid,
}

impl RigidTransform {
    /// Creates a rigid on-the-fly transform.
    #[must_use]
    pub const fn new(transform: Rigid) -> Self {
        Self { transform }
    }

    /// Creates a pure translation.
    #[must_use]
    pub const fn translation(offset: [f64; 3]) -> Self {
        Self::new(Rigid::translation(offset))
    }
}

impl FrameTransform for RigidTransform {
    fn name(&self) -> &'static str {
        "rigid"
    }

    fn apply(&self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        self.transform.apply_all(&mut timestep.positions);
        Ok(())
    }
}

/// Centres a selected set geometrically at a requested point.
#[derive(Clone, Debug)]
pub struct Center {
    atoms: Box<[usize]>,
    target: [f64; 3],
}

impl Center {
    /// Creates a geometric centre transform over explicit atom indices.
    #[must_use]
    pub fn geometric(atoms: impl Into<Box<[usize]>>, target: [f64; 3]) -> Self {
        Self {
            atoms: atoms.into(),
            target,
        }
    }
}

impl FrameTransform for Center {
    fn name(&self) -> &'static str {
        "center_geometric"
    }

    fn apply(&self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        let selected = selected_positions(&timestep.positions, &self.atoms)?;
        let Some(centre) = molframe_geom::centroid(&selected) else {
            return Ok(());
        };
        let offset = [
            self.target[0] - centre[0],
            self.target[1] - centre[1],
            self.target[2] - centre[2],
        ];
        Rigid::translation(offset).apply_all(&mut timestep.positions);
        Ok(())
    }
}

/// Fits selected mobile atoms onto fixed reference positions.
#[derive(Clone, Debug)]
pub struct Fit {
    atoms: Box<[usize]>,
    reference: Box<[[f32; 3]]>,
}

impl Fit {
    /// Creates a per-frame rigid fit with explicit correspondence.
    #[must_use]
    pub fn new(atoms: impl Into<Box<[usize]>>, reference: impl Into<Box<[[f32; 3]]>>) -> Self {
        Self {
            atoms: atoms.into(),
            reference: reference.into(),
        }
    }
}

impl FrameTransform for Fit {
    fn name(&self) -> &'static str {
        "fit"
    }

    fn apply(&self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        let mobile = selected_positions(&timestep.positions, &self.atoms)?;
        let fit = molframe_geom::superpose(&mobile, &self.reference)
            .map_err(|_| TrajectoryError::DegenerateFit)?;
        fit.transform.apply_all(&mut timestep.positions);
        Ok(())
    }
}

fn selected_positions(
    positions: &[[f32; 3]],
    atoms: &[usize],
) -> Result<Vec<[f32; 3]>, TrajectoryError> {
    atoms
        .iter()
        .map(|&index| {
            positions
                .get(index)
                .copied()
                .ok_or(TrajectoryError::SelectionOutOfRange {
                    index,
                    atoms: positions.len(),
                })
        })
        .collect()
}

/// Reader decorator applying transforms in declaration order.
#[derive(Debug)]
pub struct PipelineReader<R> {
    source: R,
    transforms: Vec<Box<dyn FrameTransform>>,
}

impl<R> PipelineReader<R> {
    /// Starts an empty transform pipeline.
    #[must_use]
    pub const fn new(source: R) -> Self {
        Self {
            source,
            transforms: Vec::new(),
        }
    }

    /// Appends one transform.
    #[must_use]
    pub fn then(mut self, transform: impl FrameTransform + 'static) -> Self {
        self.transforms.push(Box::new(transform));
        self
    }

    /// Stable transform names in application order.
    pub fn transform_names(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.transforms.iter().map(|transform| transform.name())
    }

    /// Returns the wrapped reader.
    #[must_use]
    pub fn into_inner(self) -> R {
        self.source
    }
}

impl<R: TrajectoryReader> TrajectoryReader for PipelineReader<R> {
    fn format(&self) -> &'static str {
        self.source.format()
    }

    fn n_atoms(&self) -> usize {
        self.source.n_atoms()
    }

    fn n_frames(&self) -> Option<usize> {
        self.source.n_frames()
    }

    fn units(&self) -> Units {
        self.source.units()
    }

    fn random_access(&self) -> RandomAccess {
        self.source.random_access()
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        if !self.source.read_next(timestep)? {
            return Ok(false);
        }
        for transform in &self.transforms {
            transform.apply(timestep)?;
        }
        Ok(true)
    }

    fn seek(&mut self, frame: usize) -> Result<(), TrajectoryError> {
        self.source.seek(frame)
    }
}

#[cfg(test)]
#[path = "transform_tests.rs"]
mod tests;
