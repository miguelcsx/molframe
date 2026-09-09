//! Allocation-bounded scanline flooding over compact linear cell runs.

use crate::numeric::{f64_to_f32, usize_to_f64, usize_to_u32};
use crate::workspace::empty_with_capacity;
use crate::{Cavity, SasaError};

use super::{CAVITY, EMPTY, EXTERIOR, Grid};

/// One inclusive run of cells along the contiguous X axis.
#[derive(Clone, Copy)]
struct Span {
    start: u32,
    end: u32,
}

/// One reusable scanline stack.
///
/// A checkerboard is the maximum-fragmentation case and has at most one span
/// per two cells, so the complete stack reserves four bytes per grid cell. A
/// molecular grid normally holds orders of magnitude fewer pending spans.
pub(crate) struct FloodWorkspace {
    spans: Vec<Span>,
}

impl FloodWorkspace {
    pub(crate) fn new(cell_count: usize) -> Result<Self, SasaError> {
        if cell_count > usize::try_from(u32::MAX).map_err(|_| SasaError::GridDimensionsOverflow)? {
            return Err(SasaError::GridDimensionsOverflow);
        }
        let maximum_spans = cell_count
            .checked_add(1)
            .ok_or(SasaError::GridDimensionsOverflow)?
            / 2;
        Ok(Self {
            spans: empty_with_capacity(maximum_spans)?,
        })
    }

    fn clear(&mut self) {
        self.spans.clear();
    }

    fn push(&mut self, start: usize, end: usize) {
        self.spans.push(Span {
            start: usize_to_u32(start),
            end: usize_to_u32(end),
        });
    }

    fn pop(&mut self) -> Option<(usize, usize)> {
        self.spans.pop().and_then(|span| {
            let start = usize::try_from(span.start).ok()?;
            let end = usize::try_from(span.end).ok()?;
            Some((start, end))
        })
    }
}

impl Grid {
    /// Floods all boundary-connected empty runs in `O(cells)` time.
    pub(crate) fn flood_exterior(&self, state: &mut [u8], flood: &mut FloodWorkspace) {
        let [nx, ny, nz] = self.dims;

        for z in 0..nz {
            for y in 0..ny {
                self.flood_from_seed(state, flood, self.index(0, y, z), EXTERIOR);
                self.flood_from_seed(state, flood, self.index(nx - 1, y, z), EXTERIOR);
            }
        }
        for z in 0..nz {
            for x in 1..nx - 1 {
                self.flood_from_seed(state, flood, self.index(x, 0, z), EXTERIOR);
                self.flood_from_seed(state, flood, self.index(x, ny - 1, z), EXTERIOR);
            }
        }
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                self.flood_from_seed(state, flood, self.index(x, y, 0), EXTERIOR);
                self.flood_from_seed(state, flood, self.index(x, y, nz - 1), EXTERIOR);
            }
        }
    }

    /// Labels cavities while reusing both the state array and scanline stack.
    pub(super) fn collect_cavities(
        &self,
        state: &mut [u8],
        flood: &mut FloodWorkspace,
    ) -> Result<Vec<Cavity>, SasaError> {
        let cell_volume = self.step * self.step * self.step;
        let mut cavities = empty_with_capacity(16)?;

        for index in 0..state.len() {
            if state[index] != EMPTY {
                continue;
            }

            let count = self.flood_from_seed(state, flood, index, CAVITY);
            let [x, y, z] = self.coordinates(index);
            let representative = self.centre(x, y, z);
            reserve_cavity_slot(&mut cavities)?;
            cavities.push(Cavity {
                volume: usize_to_f64(count) * cell_volume,
                representative: [
                    f64_to_f32(representative[0]),
                    f64_to_f32(representative[1]),
                    f64_to_f32(representative[2]),
                ],
                cells: count,
            });
        }
        cavities.sort_by(|left, right| right.volume.total_cmp(&left.volume));
        Ok(cavities)
    }

    /// Floods one connected component as X-axis runs, returning its cell count.
    fn flood_from_seed(
        &self,
        state: &mut [u8],
        flood: &mut FloodWorkspace,
        seed: usize,
        mark: u8,
    ) -> usize {
        if state.get(seed).copied() != Some(EMPTY) {
            return 0;
        }

        flood.clear();
        let Some((start, end)) = self.mark_run(state, seed, mark) else {
            return 0;
        };
        flood.push(start, end);
        let mut count = end - start + 1;

        while let Some((start, end)) = flood.pop() {
            count += self.expand_span(state, flood, start, end, mark);
        }
        count
    }

    /// Decodes a linear index only when a cavity representative is needed.
    fn coordinates(&self, index: usize) -> [usize; 3] {
        let nx = self.dims[0];
        let plane = nx * self.dims[1];
        let z = index / plane;
        let within_plane = index - z * plane;
        let y = within_plane / nx;
        [within_plane - y * nx, y, z]
    }

    /// Extends one seed across its complete empty X-axis run.
    fn mark_run(&self, state: &mut [u8], seed: usize, mark: u8) -> Option<(usize, usize)> {
        if state.get(seed).copied() != Some(EMPTY) {
            return None;
        }

        let nx = self.dims[0];
        let row_start = seed - seed % nx;
        let row_end = row_start + nx - 1;
        let mut start = seed;
        while start > row_start && state[start - 1] == EMPTY {
            start -= 1;
        }
        let mut end = seed;
        while end < row_end && state[end + 1] == EMPTY {
            end += 1;
        }
        state[start..=end].fill(mark);
        Some((start, end))
    }

    /// Finds and marks connected runs in the four adjacent Y/Z rows.
    fn expand_span(
        &self,
        state: &mut [u8],
        flood: &mut FloodWorkspace,
        start: usize,
        end: usize,
        mark: u8,
    ) -> usize {
        let nx = self.dims[0];
        let plane = nx * self.dims[1];
        let row_start = start - start % nx;
        let left = start - row_start;
        let right = end - row_start;
        let plane_offset = row_start % plane;
        let mut count = 0usize;

        if plane_offset >= nx {
            count += self.scan_adjacent_row(state, flood, row_start - nx, left, right, mark);
        }
        if plane_offset + nx < plane {
            count += self.scan_adjacent_row(state, flood, row_start + nx, left, right, mark);
        }
        if row_start >= plane {
            count += self.scan_adjacent_row(state, flood, row_start - plane, left, right, mark);
        }
        if row_start + plane < self.cell_count() {
            count += self.scan_adjacent_row(state, flood, row_start + plane, left, right, mark);
        }
        count
    }

    /// Scans the parent span's X interval in one adjacent row.
    fn scan_adjacent_row(
        &self,
        state: &mut [u8],
        flood: &mut FloodWorkspace,
        row_start: usize,
        left: usize,
        right: usize,
        mark: u8,
    ) -> usize {
        let mut x = left;
        let mut count = 0usize;
        while x <= right {
            let index = row_start + x;
            let Some((start, end)) = self.mark_run(state, index, mark) else {
                x += 1;
                continue;
            };
            flood.push(start, end);
            count += end - start + 1;
            x = end - row_start + 1;
        }
        count
    }
}

/// Grows the small cavity result vector through fallible, explicit requests.
fn reserve_cavity_slot(cavities: &mut Vec<Cavity>) -> Result<(), SasaError> {
    if cavities.len() < cavities.capacity() {
        return Ok(());
    }
    let additional = cavities.capacity().max(1);
    let bytes = additional
        .checked_mul(core::mem::size_of::<Cavity>())
        .ok_or(SasaError::GridDimensionsOverflow)?;
    cavities
        .try_reserve_exact(additional)
        .map_err(|_| SasaError::AllocationFailed { bytes })
}
