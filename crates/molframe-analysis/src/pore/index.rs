//! Exact atom-surface nearest-neighbour pruning for pore-grid candidates.

use super::{norm, subtract};

pub(super) const LEAF_ATOMS: usize = 16;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum IndexBuildError {
    Dimension,
    Allocation,
}

#[derive(Clone, Copy, Debug)]
struct Node {
    minimum: [f32; 3],
    maximum: [f32; 3],
    maximum_radius: f32,
    link: u32,
    count: u32,
}

impl Node {
    const fn placeholder() -> Self {
        Self {
            minimum: [f32::INFINITY; 3],
            maximum: [f32::NEG_INFINITY; 3],
            maximum_radius: 0.0,
            link: 0,
            count: 0,
        }
    }

    fn combined(left: Self, right: Self, right_index: u32) -> Self {
        Self {
            minimum: std::array::from_fn(|axis| left.minimum[axis].min(right.minimum[axis])),
            maximum: std::array::from_fn(|axis| left.maximum[axis].max(right.maximum[axis])),
            maximum_radius: left.maximum_radius.max(right.maximum_radius),
            link: right_index,
            count: 0,
        }
    }

    fn lower_bound(self, point: [f32; 3]) -> f32 {
        let delta = std::array::from_fn(|axis| {
            if point[axis] < self.minimum[axis] {
                self.minimum[axis] - point[axis]
            } else if point[axis] > self.maximum[axis] {
                point[axis] - self.maximum[axis]
            } else {
                0.0
            }
        });
        let lower_bound = norm(delta) - self.maximum_radius;
        if lower_bound.is_nan() {
            f32::NEG_INFINITY
        } else {
            // Every component is bounded using the same monotone floating-point
            // operations as the atom distance. One extra representable step
            // toward negative infinity keeps the pruning comparison conservative
            // even at a rounded atom-surface equality.
            lower_bound.next_down()
        }
    }
}

/// A balanced leaf tree whose bounds account for the largest atom radius.
#[derive(Debug)]
pub(super) struct ClearanceIndex<'a> {
    positions: &'a [[f32; 3]],
    radii: &'a [f32],
    atoms: Vec<u32>,
    nodes: Vec<Node>,
}

impl<'a> ClearanceIndex<'a> {
    pub(super) fn required_bytes(atom_count: usize) -> Option<usize> {
        let node_count = node_count(atom_count)?;
        atom_count
            .checked_mul(size_of::<u32>())?
            .checked_add(node_count.checked_mul(size_of::<Node>())?)
    }

    pub(super) fn build(
        positions: &'a [[f32; 3]],
        radii: &'a [f32],
    ) -> Result<Self, IndexBuildError> {
        let atom_count = u32::try_from(positions.len()).map_err(|_| IndexBuildError::Dimension)?;
        let expected_nodes = node_count(positions.len()).ok_or(IndexBuildError::Dimension)?;
        if u32::try_from(expected_nodes).is_err() {
            return Err(IndexBuildError::Dimension);
        }
        let mut atoms = Vec::new();
        atoms
            .try_reserve_exact(positions.len())
            .map_err(|_| IndexBuildError::Allocation)?;
        atoms.extend(0..atom_count);
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(expected_nodes)
            .map_err(|_| IndexBuildError::Allocation)?;
        let root_axis = widest_axis(positions);
        build_node(&mut atoms, 0, root_axis, positions, radii, &mut nodes)?;
        Ok(Self {
            positions,
            radii,
            atoms,
            nodes,
        })
    }

    pub(super) fn minimum_clearance(&self, point: [f32; 3], seed: u32) -> (f32, u32) {
        let seed = if usize::try_from(seed).is_ok_and(|index| index < self.positions.len()) {
            seed
        } else {
            0
        };
        let mut best_atom = seed;
        let mut best = self.atom_clearance(point, seed);
        let Some(root) = self.nodes.first().copied() else {
            return (best, best_atom);
        };
        if root.lower_bound(point) < best {
            self.search_node(0, root, point, &mut best, &mut best_atom);
        }
        (best, best_atom)
    }

    fn search_node(
        &self,
        node_index: usize,
        node: Node,
        point: [f32; 3],
        best: &mut f32,
        atom: &mut u32,
    ) {
        if node.count != 0 {
            self.search_leaf(node, point, best, atom);
            return;
        }

        let left_index = node_index + 1;
        let Ok(right_index) = usize::try_from(node.link) else {
            return;
        };
        let Some(left) = self.nodes.get(left_index).copied() else {
            return;
        };
        let Some(right) = self.nodes.get(right_index).copied() else {
            return;
        };
        let left_bound = left.lower_bound(point);
        let right_bound = right.lower_bound(point);
        if left_bound <= right_bound {
            self.search_bounded(left_index, left_bound, point, best, atom);
            self.search_bounded(right_index, right_bound, point, best, atom);
        } else {
            self.search_bounded(right_index, right_bound, point, best, atom);
            self.search_bounded(left_index, left_bound, point, best, atom);
        }
    }

    fn search_bounded(
        &self,
        node_index: usize,
        lower_bound: f32,
        point: [f32; 3],
        best: &mut f32,
        atom: &mut u32,
    ) {
        if lower_bound >= *best {
            return;
        }
        if let Some(node) = self.nodes.get(node_index).copied() {
            self.search_node(node_index, node, point, best, atom);
        }
    }

    fn search_leaf(&self, node: Node, point: [f32; 3], best: &mut f32, atom: &mut u32) {
        let (Ok(start), Ok(count)) = (usize::try_from(node.link), usize::try_from(node.count))
        else {
            return;
        };
        let Some(end) = start.checked_add(count) else {
            return;
        };
        let Some(atoms) = self.atoms.get(start..end) else {
            return;
        };
        for &candidate in atoms {
            let clearance = self.atom_clearance(point, candidate);
            if clearance < *best {
                *best = clearance;
                *atom = candidate;
            }
        }
    }

    fn atom_clearance(&self, point: [f32; 3], atom: u32) -> f32 {
        let Ok(index) = usize::try_from(atom) else {
            return f32::INFINITY;
        };
        let (Some(position), Some(radius)) = (self.positions.get(index), self.radii.get(index))
        else {
            return f32::INFINITY;
        };
        let clearance = norm(subtract(point, *position)) - radius;
        if clearance.is_nan() {
            f32::INFINITY
        } else {
            clearance
        }
    }
}

fn build_node(
    atoms: &mut [u32],
    start: usize,
    axis: usize,
    positions: &[[f32; 3]],
    radii: &[f32],
    nodes: &mut Vec<Node>,
) -> Result<u32, IndexBuildError> {
    let node_index = u32::try_from(nodes.len()).map_err(|_| IndexBuildError::Dimension)?;
    nodes.push(Node::placeholder());
    if atoms.len() <= LEAF_ATOMS {
        let node = leaf_node(atoms, start, positions, radii)?;
        nodes[usize::try_from(node_index).map_err(|_| IndexBuildError::Dimension)?] = node;
        return Ok(node_index);
    }

    let middle = atoms.len() / 2;
    atoms.select_nth_unstable_by(middle, |left, right| {
        positions[*left as usize][axis]
            .total_cmp(&positions[*right as usize][axis])
            .then_with(|| left.cmp(right))
    });
    let (left_atoms, right_atoms) = atoms.split_at_mut(middle);
    let next_axis = (axis + 1) % 3;
    let left_index = build_node(left_atoms, start, next_axis, positions, radii, nodes)?;
    if left_index != node_index + 1 {
        return Err(IndexBuildError::Dimension);
    }
    let right_start = start
        .checked_add(middle)
        .ok_or(IndexBuildError::Dimension)?;
    let right_index = build_node(right_atoms, right_start, next_axis, positions, radii, nodes)?;
    let left = nodes[usize::try_from(left_index).map_err(|_| IndexBuildError::Dimension)?];
    let right = nodes[usize::try_from(right_index).map_err(|_| IndexBuildError::Dimension)?];
    nodes[usize::try_from(node_index).map_err(|_| IndexBuildError::Dimension)?] =
        Node::combined(left, right, right_index);
    Ok(node_index)
}

fn leaf_node(
    atoms: &[u32],
    start: usize,
    positions: &[[f32; 3]],
    radii: &[f32],
) -> Result<Node, IndexBuildError> {
    let mut node = Node::placeholder();
    node.link = u32::try_from(start).map_err(|_| IndexBuildError::Dimension)?;
    node.count = u32::try_from(atoms.len()).map_err(|_| IndexBuildError::Dimension)?;
    for &atom in atoms {
        let index = usize::try_from(atom).map_err(|_| IndexBuildError::Dimension)?;
        let position = positions.get(index).ok_or(IndexBuildError::Dimension)?;
        let radius = radii.get(index).ok_or(IndexBuildError::Dimension)?;
        for (axis, &coordinate) in position.iter().enumerate() {
            node.minimum[axis] = node.minimum[axis].min(coordinate);
            node.maximum[axis] = node.maximum[axis].max(coordinate);
        }
        node.maximum_radius = node.maximum_radius.max(*radius);
    }
    Ok(node)
}

fn widest_axis(positions: &[[f32; 3]]) -> usize {
    let mut minimum = [f32::INFINITY; 3];
    let mut maximum = [f32::NEG_INFINITY; 3];
    for position in positions {
        for axis in 0..3 {
            minimum[axis] = minimum[axis].min(position[axis]);
            maximum[axis] = maximum[axis].max(position[axis]);
        }
    }
    let extent = std::array::from_fn::<_, 3, _>(|axis| maximum[axis] - minimum[axis]);
    if extent[1] > extent[0] && extent[1] >= extent[2] {
        1
    } else if extent[2] > extent[0] {
        2
    } else {
        0
    }
}

fn node_count(atom_count: usize) -> Option<usize> {
    if atom_count <= LEAF_ATOMS {
        return Some(usize::from(atom_count != 0));
    }
    let left = atom_count / 2;
    let right = atom_count - left;
    node_count(left)?
        .checked_add(node_count(right)?)?
        .checked_add(1)
}
