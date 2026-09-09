//! Compact construction of the unweighted contact graph.
//!
//! Contacts come from the shared spatial search rather than an all-pairs
//! scan, so the build is `O(N + E)` in selected sites and cutoff contacts.
//! Edges are sorted on a stable key, and components are labelled in first
//! appearance order, so the same input always yields the same graph.

use super::{NetworkBudget, NetworkError, check_memory, checked_product, checked_sum};
use pdbiox_core::ExecutionContext;
use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{
    PairQuery, PeriodicBox, SpatialError, SpatialSearchOptions, for_each_pairs_within_unsorted,
};

const INITIAL_EDGE_CAPACITY: usize = 8;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Edge {
    pub(crate) left: u32,
    pub(crate) right: u32,
}

impl Edge {
    pub(crate) fn indices(self) -> (usize, usize) {
        (expand_index(self.left), expand_index(self.right))
    }
}

#[derive(Debug)]
pub(crate) struct ContactGraph {
    pub(crate) edges: Vec<Edge>,
    pub(crate) component_of: Vec<u32>,
    pub(crate) component_sizes: Vec<u32>,
    pub(crate) maximum_degree: u32,
}

impl ContactGraph {
    pub(crate) fn build(
        positions: &[[f32; 3]],
        sites: &AtomSelection,
        selected: &[u32],
        budget: NetworkBudget,
        periodic: Option<&PeriodicBox>,
        context: &ExecutionContext,
    ) -> Result<Self, NetworkError> {
        let count = selected.len();
        let build_floor =
            checked_product(&[count, size_of::<u32>() * 3 + size_of::<u8>()], budget)?;
        check_memory(build_floor, budget)?;

        let lookup = SiteLookup::new(selected);
        let mut union = UnionFind::new(count, budget)?;
        let mut degrees = vec![0_u32; count];
        let mut edges = Vec::new();
        let mut callback_error = None;
        for_each_pairs_within_unsorted(
            &PairQuery {
                positions,
                left: sites,
                right: sites,
                cutoff: budget.contact_distance,
                options: SpatialSearchOptions::with_backend(budget.backend),
                periodic,
                context,
            },
            |pair| {
                if callback_error.is_some() {
                    return;
                }
                let (Some(left), Some(right)) =
                    (lookup.local(pair.first), lookup.local(pair.second))
                else {
                    callback_error =
                        Some(NetworkError::Spatial(SpatialError::NumericRangeExceeded));
                    return;
                };
                if left == right {
                    return;
                }
                if let Err(error) = ensure_edge_capacity(&mut edges, build_floor, budget) {
                    callback_error = Some(error);
                    return;
                }
                let (Ok(left_u32), Ok(right_u32)) = (u32::try_from(left), u32::try_from(right))
                else {
                    callback_error =
                        Some(NetworkError::Spatial(SpatialError::NumericRangeExceeded));
                    return;
                };
                edges.push(Edge {
                    left: left_u32,
                    right: right_u32,
                });
                degrees[left] += 1;
                degrees[right] += 1;
                union.join(left, right);
            },
        )?;
        if let Some(error) = callback_error {
            return Err(error);
        }
        edges.sort_unstable_by_key(|edge| (edge.left, edge.right));

        let maximum_degree = match degrees.iter().copied().max() {
            Some(maximum) => maximum,
            None => 0,
        };
        drop(degrees);
        let (component_of, component_sizes) = union.components(budget, edges.capacity())?;
        Ok(Self {
            edges,
            component_of,
            component_sizes,
            maximum_degree,
        })
    }

    pub(crate) fn site_count(&self) -> usize {
        self.component_of.len()
    }

    pub(crate) fn component_count(&self) -> usize {
        self.component_sizes.len()
    }

    pub(crate) fn owned_bytes(&self) -> usize {
        self.edges.capacity() * size_of::<Edge>()
            + self.component_of.capacity() * size_of::<u32>()
            + self.component_sizes.capacity() * size_of::<u32>()
    }
}

fn ensure_edge_capacity(
    edges: &mut Vec<Edge>,
    build_floor: usize,
    budget: NetworkBudget,
) -> Result<(), NetworkError> {
    if edges.len() < edges.capacity() {
        return Ok(());
    }
    let next = if edges.capacity() == 0 {
        INITIAL_EDGE_CAPACITY
    } else {
        edges
            .capacity()
            .checked_mul(2)
            .ok_or(NetworkError::MemoryLimit {
                required: usize::MAX,
                limit: budget.memory_limit_bytes,
            })?
    };
    let edge_bytes = checked_product(&[next, size_of::<Edge>()], budget)?;
    let required = checked_sum(&[build_floor, edge_bytes], budget)?;
    check_memory(required, budget)?;
    edges
        .try_reserve_exact(next - edges.capacity())
        .map_err(|_| NetworkError::MemoryLimit {
            required,
            limit: budget.memory_limit_bytes,
        })
}

enum SiteLookup<'a> {
    Consecutive { first: u32, len: usize },
    Sorted(&'a [u32]),
}

impl<'a> SiteLookup<'a> {
    fn new(selected: &'a [u32]) -> Self {
        let Some(&first) = selected.first() else {
            return Self::Sorted(selected);
        };
        let consecutive = selected
            .last()
            .and_then(|last| last.checked_sub(first))
            .and_then(|span| usize::try_from(span).ok())
            .and_then(|span| span.checked_add(1))
            == Some(selected.len());
        if consecutive {
            Self::Consecutive {
                first,
                len: selected.len(),
            }
        } else {
            Self::Sorted(selected)
        }
    }

    fn local(&self, atom: u32) -> Option<usize> {
        match self {
            Self::Consecutive { first, len } => atom
                .checked_sub(*first)
                .and_then(|local| usize::try_from(local).ok())
                .filter(|local| *local < *len),
            Self::Sorted(selected) => selected.binary_search(&atom).ok(),
        }
    }
}

struct UnionFind {
    parent: Vec<u32>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(count: usize, budget: NetworkBudget) -> Result<Self, NetworkError> {
        let count_u32 = u32::try_from(count).map_err(|_| NetworkError::MemoryLimit {
            required: usize::MAX,
            limit: budget.memory_limit_bytes,
        })?;
        Ok(Self {
            parent: (0..count_u32).collect(),
            rank: vec![0; count],
        })
    }

    fn root(&mut self, index: usize) -> usize {
        let mut root = index;
        while expand_index(self.parent[root]) != root {
            root = expand_index(self.parent[root]);
        }
        let mut current = index;
        while expand_index(self.parent[current]) != root {
            let next = expand_index(self.parent[current]);
            self.parent[current] = compact_index(root);
            current = next;
        }
        root
    }

    fn join(&mut self, left: usize, right: usize) {
        let mut left_root = self.root(left);
        let mut right_root = self.root(right);
        if left_root == right_root {
            return;
        }
        if self.rank[left_root] < self.rank[right_root] {
            std::mem::swap(&mut left_root, &mut right_root);
        }
        self.parent[right_root] = compact_index(left_root);
        if self.rank[left_root] == self.rank[right_root] {
            self.rank[left_root] += 1;
        }
    }

    fn components(
        self,
        budget: NetworkBudget,
        edge_capacity: usize,
    ) -> Result<(Vec<u32>, Vec<u32>), NetworkError> {
        let mut parent = self.parent;
        drop(self.rank);
        for index in 0..parent.len() {
            compress_root(&mut parent, index);
        }
        let count = parent.len();
        let component_count = parent
            .iter()
            .enumerate()
            .filter(|(index, root)| expand_index(**root) == *index)
            .count();
        let required = checked_sum(
            &[
                checked_product(&[edge_capacity, size_of::<Edge>()], budget)?,
                checked_product(&[count, size_of::<u32>() * 4], budget)?,
                checked_product(&[component_count, size_of::<u32>()], budget)?,
            ],
            budget,
        )?;
        check_memory(required, budget)?;

        let mut root_to_component = vec![u32::MAX; count];
        let mut component_sizes = vec![0_u32; component_count];
        let mut next_component = 0_u32;
        for (index, &root) in parent.iter().enumerate() {
            if expand_index(root) == index {
                root_to_component[index] = next_component;
                next_component += 1;
            }
        }
        let mut component_of = vec![0_u32; count];
        for (index, &root) in parent.iter().enumerate() {
            let component = root_to_component[expand_index(root)];
            component_of[index] = component;
            component_sizes[expand_index(component)] += 1;
        }
        Ok((component_of, component_sizes))
    }
}

fn compress_root(parent: &mut [u32], index: usize) {
    let mut root = index;
    while expand_index(parent[root]) != root {
        root = expand_index(parent[root]);
    }
    let mut current = index;
    while expand_index(parent[current]) != root {
        let next = expand_index(parent[current]);
        parent[current] = compact_index(root);
        current = next;
    }
    parent[index] = compact_index(root);
}

pub(crate) fn expand_index(index: u32) -> usize {
    match usize::try_from(index) {
        Ok(index) => index,
        Err(_) => unreachable!("a compact graph index was already validated against usize"),
    }
}

fn compact_index(index: usize) -> u32 {
    match u32::try_from(index) {
        Ok(index) => index,
        Err(_) => unreachable!("GNM site count is bounded by the public u32 atom index"),
    }
}
