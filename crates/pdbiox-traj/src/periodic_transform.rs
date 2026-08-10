//! Periodic wrapping and bond-connected molecule unwrapping.

use crate::numeric::{f32_triplet, f64_from_usize};
use crate::{FrameTransform, Timestep, TrajectoryError};
use pdbiox_spatial::PeriodicBox;
use std::collections::{BTreeSet, VecDeque};

/// Wraps atoms or explicit topology groups into the primary cell.
#[derive(Clone, Debug)]
pub struct Wrap {
    groups: Option<Vec<Box<[usize]>>>,
}

impl Wrap {
    /// Wraps every atom independently.
    #[must_use]
    pub const fn atoms() -> Self {
        Self { groups: None }
    }

    /// Wraps each explicit residue, molecule or fragment as one rigid group.
    #[must_use]
    pub fn groups(groups: impl IntoIterator<Item = impl Into<Box<[usize]>>>) -> Self {
        Self {
            groups: Some(groups.into_iter().map(Into::into).collect()),
        }
    }
}

impl FrameTransform for Wrap {
    fn name(&self) -> &'static str {
        if self.groups.is_some() {
            "wrap_groups"
        } else {
            "wrap_atoms"
        }
    }

    fn apply(&self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        let periodic = periodic(timestep)?;
        match &self.groups {
            None => {
                for position in &mut timestep.positions {
                    *position = periodic.wrap(*position);
                }
            }
            Some(groups) => wrap_groups(&mut timestep.positions, groups, periodic)?,
        }
        Ok(())
    }
}

fn wrap_groups(
    positions: &mut [[f32; 3]],
    groups: &[Box<[usize]>],
    periodic: PeriodicBox,
) -> Result<(), TrajectoryError> {
    let mut seen = BTreeSet::new();
    for group in groups {
        let mut centre = [0.0_f64; 3];
        for &index in group {
            let Some(position) = positions.get(index) else {
                return Err(out_of_range(index, positions.len()));
            };
            if !seen.insert(index) {
                return Err(TrajectoryError::OverlappingGroups { index });
            }
            for axis in 0..3 {
                centre[axis] += f64::from(position[axis]);
            }
        }
        if group.is_empty() {
            continue;
        }
        let divisor =
            f64_from_usize(group.len()).ok_or(TrajectoryError::UnrepresentableCoordinate)?;
        for value in &mut centre {
            *value /= divisor;
        }
        let centre = f32_triplet(centre).ok_or(TrajectoryError::UnrepresentableCoordinate)?;
        let fractional = periodic.fractional(centre);
        let shift = periodic.cartesian(fractional.map(|value| -value.floor()));
        for &index in group {
            for (coordinate, offset) in positions[index].iter_mut().zip(shift) {
                *coordinate += offset;
            }
        }
    }
    Ok(())
}

/// Makes bond-connected molecules whole across periodic boundaries.
#[derive(Clone, Debug)]
pub struct Unwrap {
    bonds: Box<[(usize, usize)]>,
}

impl Unwrap {
    /// Creates a molecule-unwrapping transform from topology bonds.
    #[must_use]
    pub fn molecules(bonds: impl Into<Box<[(usize, usize)]>>) -> Self {
        Self {
            bonds: bonds.into(),
        }
    }
}

impl FrameTransform for Unwrap {
    fn name(&self) -> &'static str {
        "make_molecules_whole"
    }

    fn apply(&self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        let periodic = periodic(timestep)?;
        let atom_count = timestep.positions.len();
        let mut adjacency = vec![Vec::new(); atom_count];
        for &(left, right) in &self.bonds {
            if left >= atom_count {
                return Err(out_of_range(left, atom_count));
            }
            if right >= atom_count {
                return Err(out_of_range(right, atom_count));
            }
            adjacency[left].push(right);
            adjacency[right].push(left);
        }
        for neighbors in &mut adjacency {
            neighbors.sort_unstable();
            neighbors.dedup();
        }
        let original = timestep.positions.clone();
        let mut visited = vec![false; atom_count];
        for root in 0..atom_count {
            if visited[root] || adjacency[root].is_empty() {
                continue;
            }
            visited[root] = true;
            let mut queue = VecDeque::from([root]);
            while let Some(parent) = queue.pop_front() {
                for &child in &adjacency[parent] {
                    if visited[child] {
                        continue;
                    }
                    let displacement = periodic.displacement(original[parent], original[child]);
                    timestep.positions[child] = [
                        timestep.positions[parent][0] + displacement[0],
                        timestep.positions[parent][1] + displacement[1],
                        timestep.positions[parent][2] + displacement[2],
                    ];
                    visited[child] = true;
                    queue.push_back(child);
                }
            }
        }
        Ok(())
    }
}

fn periodic(timestep: &Timestep) -> Result<PeriodicBox, TrajectoryError> {
    timestep
        .cell
        .ok_or(TrajectoryError::MissingCell)
        .and_then(|cell| PeriodicBox::from_cell(cell).map_err(Into::into))
}

const fn out_of_range(index: usize, atoms: usize) -> TrajectoryError {
    TrajectoryError::SelectionOutOfRange { index, atoms }
}

#[cfg(test)]
#[path = "periodic_transform_tests.rs"]
mod tests;
