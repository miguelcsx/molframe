//! Dense micro-kernel and restarted matrix-free solve.

use super::{GaussianNetworkModel, GnmError, GnmOptions, budget, canonical_mode};
use crate::network::{self, ContactGraph, expand_index};
use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::matrix_free::LinOp;
use faer::matrix_free::eigen::{
    PartialEigenParams, partial_self_adjoint_eigen, partial_self_adjoint_eigen_scratch,
};
use faer::reborrow::ReborrowMut;
use faer::{Col, Mat, MatMut, MatRef, Par};
use nalgebra::{DMatrix, SymmetricEigen};

const DENSE_SITE_MAXIMUM: usize = 64;
const DENSE_MATRIX_COUNT: usize = 2;
const SOLVER_RESTARTS: usize = 1_000;

/// One Gaussian mode exists per site, less one rigid translation per connected
/// component.
fn available_modes(graph: &ContactGraph) -> usize {
    graph.site_count() - graph.component_count()
}

pub(super) fn solve(
    graph: &ContactGraph,
    sites: Vec<u32>,
    options: GnmOptions,
) -> Result<GaussianNetworkModel, GnmError> {
    let available = available_modes(graph);
    if available < options.mode_count {
        return Err(GnmError::InsufficientModes {
            requested: options.mode_count,
            available,
        });
    }
    if graph.site_count() <= DENSE_SITE_MAXIMUM {
        dense_modes(graph, sites, options)
    } else {
        sparse_modes(graph, sites, options)
    }
}

fn dense_modes(
    graph: &ContactGraph,
    sites: Vec<u32>,
    options: GnmOptions,
) -> Result<GaussianNetworkModel, GnmError> {
    let count = graph.site_count();
    let matrix_bytes = checked_product(
        &[count, count, DENSE_MATRIX_COUNT, size_of::<f64>()],
        options,
    )?;
    let output_bytes = checked_product(&[count, options.mode_count, size_of::<f64>()], options)?;
    let site_bytes = checked_product(&[sites.capacity(), size_of::<u32>()], options)?;
    let required = checked_sum(
        &[graph.owned_bytes(), site_bytes, matrix_bytes, output_bytes],
        options,
    )?;
    check_memory(required, options)?;

    let mut kirchhoff = DMatrix::zeros(count, count);
    for edge in &graph.edges {
        let (left, right) = edge.indices();
        kirchhoff[(left, right)] = -1.0;
        kirchhoff[(right, left)] = -1.0;
        kirchhoff[(left, left)] += 1.0;
        kirchhoff[(right, right)] += 1.0;
    }
    modes_from_matrix(kirchhoff, sites, options)
}

fn sparse_modes(
    graph: &ContactGraph,
    sites: Vec<u32>,
    options: GnmOptions,
) -> Result<GaussianNetworkModel, GnmError> {
    let available = available_modes(graph);
    let mut solve_count = options.mode_count;
    loop {
        let candidates = solve_sparse_candidates(graph, &sites, solve_count, options)?;
        let numerical_zero_count = candidates
            .iter()
            .take_while(|candidate| candidate.value.abs() <= options.zero_mode_tolerance)
            .count();
        let non_zero_count = candidates.len() - numerical_zero_count;
        if non_zero_count >= options.mode_count {
            let chosen =
                &candidates[numerical_zero_count..numerical_zero_count + options.mode_count];
            return Ok(GaussianNetworkModel {
                sites,
                eigenvalues: chosen.iter().map(|candidate| candidate.value).collect(),
                modes: chosen
                    .iter()
                    .map(|candidate| candidate.mode.clone())
                    .collect(),
                zero_modes: graph.component_count() + numerical_zero_count,
            });
        }
        if solve_count == available {
            return Err(GnmError::InsufficientModes {
                requested: options.mode_count,
                available: non_zero_count,
            });
        }
        let missing = options.mode_count - non_zero_count;
        let doubled = solve_count.checked_mul(2).ok_or(GnmError::MemoryLimit {
            required: usize::MAX,
            limit: options.memory_limit_bytes,
        })?;
        let filled = solve_count
            .checked_add(missing)
            .ok_or(GnmError::MemoryLimit {
                required: usize::MAX,
                limit: options.memory_limit_bytes,
            })?;
        solve_count = available.min(doubled.max(filled));
    }
}

#[derive(Debug)]
struct Candidate {
    value: f64,
    shifted_value: f64,
    mode: Vec<f64>,
}

fn solve_sparse_candidates(
    graph: &ContactGraph,
    sites: &[u32],
    solve_count: usize,
    options: GnmOptions,
) -> Result<Vec<Candidate>, GnmError> {
    let shift = 2.0 * f64::from(graph.maximum_degree) + 1.0;
    let params = PartialEigenParams {
        max_restarts: SOLVER_RESTARTS,
        ..PartialEigenParams::default()
    };
    // The solver's private Rayon path would bypass ExecutionContext and can
    // oversubscribe nested calls, so this kernel remains sequential until it
    // can submit bounded work to the shared executor.
    let par = Par::Seq;
    let mut candidates = Vec::with_capacity(solve_count);
    let maximum_batch = (graph.site_count() - 1) / 2;
    while candidates.len() < solve_count {
        let remaining = solve_count - candidates.len();
        let batch_count = if graph.component_count() == 1 {
            remaining.min(maximum_batch)
        } else {
            1
        };
        let operator = ShiftedLaplacian {
            graph,
            shift,
            deflated: &candidates,
        };
        let scratch =
            partial_self_adjoint_eigen_scratch::<f64>(&operator, batch_count, par, params);
        check_sparse_memory(graph, sites, solve_count, scratch, options)?;

        let tolerance = f64::EPSILON * 256.0 * shift;
        let initial = deterministic_initial_vector(graph, &candidates);
        if let Some(candidate) = exact_candidate(&operator, &initial, tolerance) {
            candidates.push(candidate);
            continue;
        }
        let mut eigenvectors = Mat::<f64>::zeros(graph.site_count(), batch_count);
        let mut shifted_values = vec![0.0; batch_count];
        let mut memory = MemBuffer::new(scratch);
        let info = partial_self_adjoint_eigen(
            eigenvectors.as_mut(),
            &mut shifted_values,
            &operator,
            initial.as_ref(),
            tolerance,
            par,
            MemStack::new(&mut memory),
            params,
        );
        if info.n_converged_eigen == 0 {
            return Err(GnmError::Convergence {
                requested: solve_count,
                converged: candidates.len(),
            });
        }
        for (column, _) in shifted_values
            .into_iter()
            .take(info.n_converged_eigen)
            .enumerate()
        {
            let mode = canonical_mode(eigenvectors.col(column).iter().copied().collect());
            let value = laplacian_rayleigh(graph, &mode);
            candidates.push(Candidate {
                value,
                shifted_value: shift - value,
                mode,
            });
        }
    }
    candidates.sort_by(|left, right| left.value.total_cmp(&right.value));
    Ok(candidates)
}

fn exact_candidate(op: &ShiftedLaplacian<'_>, initial: &Col<f64>, tol: f64) -> Option<Candidate> {
    let mut applied = Col::zeros(op.graph.site_count());
    let mut memory = MemBuffer::new(op.apply_scratch(1, Par::Seq));
    let stack = MemStack::new(&mut memory);
    op.apply(applied.as_mat_mut(), initial.as_mat(), Par::Seq, stack);
    let norm: f64 = initial.iter().map(|value| value * value).sum();
    let shifted = initial
        .iter()
        .zip(applied.iter())
        .map(|(x, y)| x * y)
        .sum::<f64>()
        / norm;
    let residual: f64 = initial
        .iter()
        .zip(applied.iter())
        .map(|(x, y)| (y - shifted * x).powi(2))
        .sum::<f64>()
        / norm;
    (residual.sqrt() <= tol).then(|| {
        let mode = canonical_mode(initial.iter().map(|value| value / norm.sqrt()).collect());
        let value = laplacian_rayleigh(op.graph, &mode);
        Candidate {
            value,
            shifted_value: op.shift - value,
            mode,
        }
    })
}

fn laplacian_rayleigh(graph: &ContactGraph, mode: &[f64]) -> f64 {
    let squared_norm: f64 = mode.iter().map(|value| value * value).sum();
    let energy: f64 = graph
        .edges
        .iter()
        .map(|edge| {
            let (left, right) = edge.indices();
            let difference = mode[left] - mode[right];
            difference * difference
        })
        .sum();
    energy / squared_norm
}

fn check_sparse_memory(
    graph: &ContactGraph,
    sites: &[u32],
    solve_count: usize,
    scratch: StackReq,
    options: GnmOptions,
) -> Result<(), GnmError> {
    let vector_bytes = checked_product(
        &[graph.site_count(), solve_count, size_of::<f64>()],
        options,
    )?;
    let doubled_vectors = checked_product(&[vector_bytes, 2], options)?;
    let site_bytes = checked_product(&[sites.len(), size_of::<u32>()], options)?;
    let initial_bytes = checked_product(&[graph.site_count(), size_of::<f64>()], options)?;
    let value_bytes = checked_product(&[solve_count, size_of::<f64>()], options)?;
    let candidate_bytes = checked_product(&[solve_count, size_of::<Candidate>()], options)?;
    let required = checked_sum(
        &[
            graph.owned_bytes(),
            site_bytes,
            doubled_vectors,
            initial_bytes,
            value_bytes,
            candidate_bytes,
            scratch.unaligned_bytes_required(),
        ],
        options,
    )?;
    check_memory(required, options)
}

fn deterministic_initial_vector(graph: &ContactGraph, deflated: &[Candidate]) -> Col<f64> {
    let mut vector = Col::from_fn(graph.site_count(), |index| {
        let Ok(index) = u32::try_from(index) else {
            unreachable!("GNM site count is bounded by the public atom index")
        };
        let Ok(deflated_count) = u32::try_from(deflated.len()) else {
            unreachable!("GNM mode count cannot exceed its u32 site count")
        };
        let mut value = u64::from(index)
            .wrapping_add(0x9e37_79b9_7f4a_7c15)
            .wrapping_add(u64::from(deflated_count).wrapping_mul(0xd1b5_4a32_d192_ed03));
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^= value >> 31;
        let Ok(upper) = u32::try_from(value >> 32) else {
            unreachable!("the shifted hash occupies exactly 32 bits")
        };
        (f64::from(upper) / f64::from(u32::MAX)) - 0.5
    });
    let mut sums = vec![0.0; graph.component_count()];
    for (index, &component) in graph.component_of.iter().enumerate() {
        sums[expand_index(component)] += vector[index];
    }
    for (index, &component) in graph.component_of.iter().enumerate() {
        let component = expand_index(component);
        vector[index] -= sums[component] / f64::from(graph.component_sizes[component]);
    }
    for candidate in deflated {
        let projection: f64 = candidate
            .mode
            .iter()
            .enumerate()
            .map(|(index, value)| value * vector[index])
            .sum();
        for (value, basis) in vector.iter_mut().zip(&candidate.mode) {
            *value -= projection * basis;
        }
    }
    vector
}

#[derive(Debug)]
struct ShiftedLaplacian<'a> {
    graph: &'a ContactGraph,
    shift: f64,
    deflated: &'a [Candidate],
}

impl LinOp<f64> for ShiftedLaplacian<'_> {
    fn apply_scratch(&self, _rhs_ncols: usize, _par: Par) -> StackReq {
        StackReq::new::<f64>(self.graph.component_count())
    }

    fn nrows(&self) -> usize {
        self.graph.site_count()
    }

    fn ncols(&self) -> usize {
        self.graph.site_count()
    }

    fn apply(
        &self,
        mut out: MatMut<'_, f64>,
        rhs: MatRef<'_, f64>,
        _par: Par,
        stack: &mut MemStack,
    ) {
        for column in 0..rhs.ncols() {
            let out_column = out.rb_mut().col_mut(column).try_as_col_major_mut();
            let rhs_column = rhs.col(column).try_as_col_major();
            match (out_column, rhs_column) {
                (Some(out_column), Some(rhs_column)) => {
                    self.apply_contiguous(out_column.as_slice_mut(), rhs_column.as_slice(), stack);
                }
                _ => self.apply_strided(out.rb_mut(), rhs, column, stack),
            }
        }
    }

    fn conj_apply(
        &self,
        out: MatMut<'_, f64>,
        rhs: MatRef<'_, f64>,
        par: Par,
        stack: &mut MemStack,
    ) {
        self.apply(out, rhs, par, stack);
    }
}

impl ShiftedLaplacian<'_> {
    fn apply_contiguous(&self, out: &mut [f64], rhs: &[f64], stack: &mut MemStack) {
        let (mut component_sums, _) = stack.make_with(self.graph.component_count(), |_| 0.0_f64);
        for ((out, rhs), &component) in out.iter_mut().zip(rhs).zip(&self.graph.component_of) {
            *out = self.shift * *rhs;
            component_sums[expand_index(component)] += rhs;
        }
        for (out, &component) in out.iter_mut().zip(&self.graph.component_of) {
            let component = expand_index(component);
            let size = self.graph.component_sizes[component];
            *out -= self.shift * component_sums[component] / f64::from(size);
        }
        for edge in &self.graph.edges {
            let (left, right) = edge.indices();
            let difference = rhs[left] - rhs[right];
            out[left] -= difference;
            out[right] += difference;
        }
        for candidate in self.deflated {
            let projection: f64 = candidate
                .mode
                .iter()
                .zip(rhs)
                .map(|(basis, rhs)| basis * rhs)
                .sum();
            for (out, basis) in out.iter_mut().zip(&candidate.mode) {
                *out -= candidate.shifted_value * projection * basis;
            }
        }
    }

    fn apply_strided(
        &self,
        mut out: MatMut<'_, f64>,
        rhs: MatRef<'_, f64>,
        column: usize,
        stack: &mut MemStack,
    ) {
        let (mut component_sums, _) = stack.make_with(self.graph.component_count(), |_| 0.0_f64);
        for (row, &component) in self.graph.component_of.iter().enumerate() {
            out[(row, column)] = self.shift * rhs[(row, column)];
            component_sums[expand_index(component)] += rhs[(row, column)];
        }
        for (row, &component) in self.graph.component_of.iter().enumerate() {
            let component = expand_index(component);
            let size = self.graph.component_sizes[component];
            out[(row, column)] -= self.shift * component_sums[component] / f64::from(size);
        }
        for edge in &self.graph.edges {
            let (left, right) = edge.indices();
            let difference = rhs[(left, column)] - rhs[(right, column)];
            out[(left, column)] -= difference;
            out[(right, column)] += difference;
        }
        for candidate in self.deflated {
            let projection: f64 = candidate
                .mode
                .iter()
                .enumerate()
                .map(|(row, value)| value * rhs[(row, column)])
                .sum();
            for (row, value) in candidate.mode.iter().enumerate() {
                out[(row, column)] -= candidate.shifted_value * projection * value;
            }
        }
    }
}

pub(super) fn modes_from_matrix(
    kirchhoff: DMatrix<f64>,
    sites: Vec<u32>,
    options: GnmOptions,
) -> Result<GaussianNetworkModel, GnmError> {
    let eigen = SymmetricEigen::new(kirchhoff);
    let mut indexed: Vec<(f64, usize)> = eigen
        .eigenvalues
        .iter()
        .copied()
        .enumerate()
        .map(|(index, value)| (value, index))
        .collect();
    indexed.sort_by(|left, right| left.0.total_cmp(&right.0));
    let zero_modes = indexed
        .iter()
        .take_while(|(value, _)| value.abs() <= options.zero_mode_tolerance)
        .count();
    let available = indexed.len() - zero_modes;
    if available < options.mode_count {
        return Err(GnmError::InsufficientModes {
            requested: options.mode_count,
            available,
        });
    }
    let chosen = &indexed[zero_modes..zero_modes + options.mode_count];
    let eigenvalues = chosen.iter().map(|(value, _)| *value).collect();
    let modes = chosen
        .iter()
        .map(|(_, column)| {
            canonical_mode(eigen.eigenvectors.column(*column).iter().copied().collect())
        })
        .collect();
    Ok(GaussianNetworkModel {
        sites,
        eigenvalues,
        modes,
        zero_modes,
    })
}

pub(super) fn checked_product(factors: &[usize], options: GnmOptions) -> Result<usize, GnmError> {
    network::checked_product(factors, budget(options)).map_err(GnmError::from)
}

pub(super) fn checked_sum(values: &[usize], options: GnmOptions) -> Result<usize, GnmError> {
    network::checked_sum(values, budget(options)).map_err(GnmError::from)
}

pub(super) fn check_memory(required: usize, options: GnmOptions) -> Result<(), GnmError> {
    network::check_memory(required, budget(options)).map_err(GnmError::from)
}
