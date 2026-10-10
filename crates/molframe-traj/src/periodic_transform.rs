//! Periodic wrapping and bond-connected molecule un-wrapping.

use crate::numeric::{f32_from_f64, f32_triplet, f64_from_usize};
use crate::{FrameTransform, Timestep, TrajectoryError};
use molframe_spatial::PeriodicBox;
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

    fn apply(&mut self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
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
///
/// This is a per-frame operation: each molecule is reassembled around one of
/// its own atoms so no bond spans the box, but nothing relates one frame to the
/// next. A molecule diffusing across a boundary still jumps by a box vector
/// between frames. Use [`TemporalUnwrap`] when continuity over time matters
/// (diffusion, mean-squared displacement, centre-of-mass tracking).
#[derive(Clone, Debug)]
pub struct Unwrap {
    bonds: Box<[(usize, usize)]>,
}

impl Unwrap {
    /// Creates a make-molecules-whole transform from topology bonds.
    ///
    /// Same as [`Unwrap::make_molecules_whole`]; the name is kept for
    /// compatibility.
    #[must_use]
    pub fn molecules(bonds: impl Into<Box<[(usize, usize)]>>) -> Self {
        Self {
            bonds: bonds.into(),
        }
    }

    /// Creates a make-molecules-whole transform from topology bonds.
    ///
    /// The result is whole within each frame only; it is not continuous across
    /// frames.
    #[must_use]
    pub fn make_molecules_whole(bonds: impl Into<Box<[(usize, usize)]>>) -> Self {
        Self::molecules(bonds)
    }
}

impl FrameTransform for Unwrap {
    fn name(&self) -> &'static str {
        "make_molecules_whole"
    }

    fn apply(&mut self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        let periodic = periodic(timestep)?;
        make_whole(&mut timestep.positions, &self.bonds, periodic)
    }
}

fn make_whole(
    positions: &mut [[f32; 3]],
    bonds: &[(usize, usize)],
    periodic: PeriodicBox,
) -> Result<(), TrajectoryError> {
    let atom_count = positions.len();
    let mut adjacency = vec![Vec::new(); atom_count];
    for &(left, right) in bonds {
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
    let original = positions.to_vec();
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
                positions[child] = [
                    positions[parent][0] + displacement[0],
                    positions[parent][1] + displacement[1],
                    positions[parent][2] + displacement[2],
                ];
                visited[child] = true;
                queue.push_back(child);
            }
        }
    }
    Ok(())
}

/// Continuous-in-time unwrapping of atoms or bonded molecules.
///
/// Each unit (an atom, or a bond-connected molecule) is moved so that its
/// position or geometric centre changes between consecutive frames by the
/// minimum-image displacement of its wrapped motion. A particle that leaves
/// through one face therefore continues past the box edge instead of
/// reappearing on the other side. Molecules are first made whole within the
/// frame, so internal geometry is preserved.
///
/// The transform is stateful: it must see frames in order, starting from the
/// first frame to be analysed, and it assumes no unit moves more than half a
/// box between consecutive frames. The first frame it sees is left as is.
/// Construct a fresh one per pass over a trajectory.
#[derive(Clone, Debug)]
pub struct TemporalUnwrap {
    bonds: Box<[(usize, usize)]>,
    state: Option<UnwrapState>,
}

/// Per-unit centres from the previous frame.
#[derive(Clone, Debug)]
struct UnwrapState {
    /// Unit membership, derived from the bonds on the first frame.
    units: Vec<Box<[usize]>>,
    /// Centre of each unit before temporal correction in the previous frame.
    previous_raw: Vec<[f64; 3]>,
    /// Continuous centre of each unit in the previous frame.
    previous_continuous: Vec<[f64; 3]>,
}

impl TemporalUnwrap {
    /// Tracks every atom independently.
    #[must_use]
    pub fn atoms() -> Self {
        Self {
            bonds: Box::new([]),
            state: None,
        }
    }

    /// Tracks each bond-connected molecule by its geometric centre (the mean of
    /// its atom positions; masses are not available here), keeping it whole in
    /// every frame.
    ///
    /// Atoms without bonds are tracked individually.
    #[must_use]
    pub fn molecules(bonds: impl Into<Box<[(usize, usize)]>>) -> Self {
        Self {
            bonds: bonds.into(),
            state: None,
        }
    }
}

impl FrameTransform for TemporalUnwrap {
    fn name(&self) -> &'static str {
        "temporal_unwrap"
    }

    fn apply(&mut self, timestep: &mut Timestep) -> Result<(), TrajectoryError> {
        let periodic = periodic(timestep)?;
        let count = timestep.positions.len();
        if !self.bonds.is_empty() {
            make_whole(&mut timestep.positions, &self.bonds, periodic)?;
        }
        let state = match self.state.take() {
            Some(state) if state_covers(&state, count) => state,
            _ => UnwrapState {
                units: units(count, &self.bonds)?,
                previous_raw: Vec::new(),
                previous_continuous: Vec::new(),
            },
        };
        let mut raw = Vec::with_capacity(state.units.len());
        let mut continuous = Vec::with_capacity(state.units.len());
        for (position, unit) in state.units.iter().enumerate() {
            let centre = centroid(&timestep.positions, unit)?;
            let target = match (
                state.previous_raw.get(position),
                state.previous_continuous.get(position),
            ) {
                (Some(before), Some(before_continuous)) => {
                    let step = periodic.displacement_f64(*before, centre);
                    [
                        before_continuous[0] + step[0],
                        before_continuous[1] + step[1],
                        before_continuous[2] + step[2],
                    ]
                }
                _ => centre,
            };
            for &index in unit {
                for axis in 0..3 {
                    let moved =
                        f64::from(timestep.positions[index][axis]) + target[axis] - centre[axis];
                    timestep.positions[index][axis] =
                        f32_from_f64(moved).ok_or(TrajectoryError::UnrepresentableCoordinate)?;
                }
            }
            raw.push(centre);
            continuous.push(target);
        }
        self.state = Some(UnwrapState {
            units: state.units,
            previous_raw: raw,
            previous_continuous: continuous,
        });
        Ok(())
    }
}

fn state_covers(state: &UnwrapState, atoms: usize) -> bool {
    state.units.iter().map(|unit| unit.len()).sum::<usize>() == atoms
}

fn find(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

/// Connected components of the bond graph, singletons included, by lowest atom.
fn units(count: usize, bonds: &[(usize, usize)]) -> Result<Vec<Box<[usize]>>, TrajectoryError> {
    let mut parent: Vec<usize> = (0..count).collect();
    for &(left, right) in bonds {
        if left >= count {
            return Err(out_of_range(left, count));
        }
        if right >= count {
            return Err(out_of_range(right, count));
        }
        let (a, b) = (find(&mut parent, left), find(&mut parent, right));
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> =
        std::collections::BTreeMap::new();
    for atom in 0..count {
        let root = find(&mut parent, atom);
        groups.entry(root).or_default().push(atom);
    }
    Ok(groups.into_values().map(Vec::into_boxed_slice).collect())
}

fn centroid(positions: &[[f32; 3]], unit: &[usize]) -> Result<[f64; 3], TrajectoryError> {
    let mut sum = [0.0_f64; 3];
    for &index in unit {
        for axis in 0..3 {
            sum[axis] += f64::from(positions[index][axis]);
        }
    }
    let divisor = f64_from_usize(unit.len()).ok_or(TrajectoryError::UnrepresentableCoordinate)?;
    Ok(sum.map(|value| value / divisor))
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
