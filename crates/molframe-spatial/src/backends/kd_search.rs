//! Iterative k-d traversal, bounded nearest-neighbour heap and periodic images.
//!
//! Traversal keeps an explicit stack of frames carrying the per-axis offset
//! from the query to the node's cell (Arya and Mount's incremental distance),
//! so the cell distance is three multiplies and no tree node is revisited.

use super::heap::KnnHeap;
use super::{KdTree, LEAF_MAX};
use crate::SpatialError;
use crate::brute::finite;

/// Relative slack applied to pruning so that the incremental cell distance,
/// which rounds differently from the exact point distance, never discards a
/// point whose exact squared distance is within the bound.
const PRUNE_SLACK: f32 = 1.0 + 1.0e-5;

/// Image shifts visited for an unbounded periodic nearest-neighbour query.
pub(super) const NEAREST_IMAGE_LIMITS: [i64; 3] = [2; 3];

#[derive(Clone, Copy, Debug)]
pub(super) struct Frame {
    node: u32,
    lo: u32,
    hi: u32,
    off: [f32; 3],
    rd: f32,
}

/// Receives leaf scans and tells the traversal how far it still needs to look.
pub(super) trait Sink {
    /// Squared distance beyond which no further point can be accepted.
    fn bound(&self) -> f32;
    /// Offers one leaf bucket: atom ids and their squared distances.
    fn leaf(&mut self, ids: &[u32], squared: &[f32]);
}

/// Reusable per-worker traversal storage; queries never allocate after warm-up.
#[derive(Debug, Default)]
pub(super) struct Scratch {
    stack: Vec<Frame>,
    candidates: Vec<u32>,
    exact: Vec<(f32, u32)>,
}

impl KdTree<'_> {
    /// Squared distance from `point` to the tree's bounding box, per axis.
    fn box_offsets(&self, point: [f32; 3]) -> Option<([f32; 3], f32)> {
        let (low, high) = self.bounds?;
        let mut off = [0.0_f32; 3];
        for axis in 0..3 {
            if point[axis] < low[axis] {
                off[axis] = point[axis] - low[axis];
            } else if point[axis] > high[axis] {
                off[axis] = point[axis] - high[axis];
            }
        }
        Some((off, off[0] * off[0] + off[1] * off[1] + off[2] * off[2]))
    }

    /// Visits every leaf bucket that can hold a point within the sink bound.
    pub(super) fn traverse<S: Sink>(&self, point: [f32; 3], stack: &mut Vec<Frame>, sink: &mut S) {
        let Some((off, rd)) = self.box_offsets(point) else {
            return;
        };
        let Ok(count) = u32::try_from(self.ids.len()) else {
            return;
        };
        let internal = self.split_value.len();
        let mut buffer = [0.0_f32; LEAF_MAX];
        stack.clear();
        stack.push(Frame {
            node: 0,
            lo: 0,
            hi: count,
            off,
            rd,
        });
        while let Some(mut frame) = stack.pop() {
            loop {
                if frame.rd > sink.bound() * PRUNE_SLACK {
                    break;
                }
                let node = frame.node as usize;
                if node >= internal {
                    self.scan_leaf(
                        point,
                        frame.lo as usize,
                        frame.hi as usize,
                        &mut buffer,
                        sink,
                    );
                    break;
                }
                let axis = usize::from(self.split_axis[node]);
                let delta = point[axis] - self.split_value[node];
                let middle = frame.lo + (frame.hi - frame.lo) / 2;
                let left = Frame {
                    node: frame.node * 2 + 1,
                    lo: frame.lo,
                    hi: middle,
                    ..frame
                };
                let right = Frame {
                    node: frame.node * 2 + 2,
                    lo: middle,
                    hi: frame.hi,
                    ..frame
                };
                let (near, mut far) = if delta < 0.0 {
                    (left, right)
                } else {
                    (right, left)
                };
                far.off[axis] = delta;
                far.rd =
                    far.off[0] * far.off[0] + far.off[1] * far.off[1] + far.off[2] * far.off[2];
                if far.rd <= sink.bound() * PRUNE_SLACK {
                    stack.push(far);
                }
                frame = near;
            }
        }
    }

    #[inline]
    fn scan_leaf<S: Sink>(
        &self,
        point: [f32; 3],
        lo: usize,
        hi: usize,
        buffer: &mut [f32; LEAF_MAX],
        sink: &mut S,
    ) {
        let (Some(xs), Some(ys), Some(zs), Some(ids)) = (
            self.xs.get(lo..hi),
            self.ys.get(lo..hi),
            self.zs.get(lo..hi),
            self.ids.get(lo..hi),
        ) else {
            return;
        };
        let squared = &mut buffer[..ids.len().min(LEAF_MAX)];
        for (((out, &x), &y), &z) in squared.iter_mut().zip(xs).zip(ys).zip(zs) {
            let dx = point[0] - x;
            let dy = point[1] - y;
            let dz = point[2] - z;
            *out = dx * dx + dy * dy + dz * dz;
        }
        sink.leaf(ids, squared);
    }

    /// Image translations ordered nearest first, as exact Cartesian offsets.
    pub(super) fn image_translations(
        &self,
        limits: [i64; 3],
        out: &mut Vec<[f64; 3]>,
    ) -> Result<(), SpatialError> {
        out.clear();
        let Some(periodic) = self.periodic else {
            return Ok(());
        };
        super::for_each_shift(limits, |shift| {
            let origin = periodic.translated([0.0; 3], shift)?;
            out.push(origin.map(f64::from));
            Ok(())
        })?;
        out.sort_by(|a, b| {
            let norm = |t: &[f64; 3]| t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
            norm(a).total_cmp(&norm(b))
        });
        Ok(())
    }

    /// Runs `sink` over the query's tree images (one image when non-periodic).
    pub(super) fn search_images<S: Sink>(
        &self,
        point: [f32; 3],
        images: &[[f64; 3]],
        stack: &mut Vec<Frame>,
        sink: &mut S,
    ) {
        let Some(periodic) = self.periodic else {
            self.traverse(point, stack, sink);
            return;
        };
        let wrapped = periodic.wrap(point);
        #[allow(clippy::cast_possible_truncation)]
        for shift in images {
            let image = [
                (f64::from(wrapped[0]) + shift[0]) as f32,
                (f64::from(wrapped[1]) + shift[1]) as f32,
                (f64::from(wrapped[2]) + shift[2]) as f32,
            ];
            self.traverse(image, stack, sink);
        }
    }
}

/// Collects every point within a bound, optionally excluding one atom.
pub(super) struct RadiusSink<'e, F: FnMut(u32, f32)> {
    pub(super) bound: f32,
    pub(super) exclude: Option<u32>,
    pub(super) emit: &'e mut F,
}

impl<F: FnMut(u32, f32)> Sink for RadiusSink<'_, F> {
    #[inline]
    fn bound(&self) -> f32 {
        self.bound
    }

    #[inline]
    fn leaf(&mut self, ids: &[u32], squared: &[f32]) {
        for (&id, &d2) in ids.iter().zip(squared) {
            if d2 <= self.bound && Some(id) != self.exclude {
                (self.emit)(id, d2);
            }
        }
    }
}

/// Feeds leaf scans into a [`KnnHeap`], optionally skipping one atom.
pub(super) struct KnnSink<'h> {
    pub(super) heap: &'h mut KnnHeap,
    pub(super) exclude: Option<u32>,
}

impl Sink for KnnSink<'_> {
    #[inline]
    fn bound(&self) -> f32 {
        self.heap.current_bound()
    }

    #[inline]
    fn leaf(&mut self, ids: &[u32], squared: &[f32]) {
        for (&id, &d2) in ids.iter().zip(squared) {
            if d2 <= self.heap.current_bound() && Some(id) != self.exclude {
                self.heap.push(d2, id);
            }
        }
    }
}

/// Collects candidate atom ids within an inflated bound.
struct CandidateSink<'v> {
    bound: f32,
    exclude: Option<u32>,
    ids: &'v mut Vec<u32>,
}

impl Sink for CandidateSink<'_> {
    #[inline]
    fn bound(&self) -> f32 {
        self.bound
    }

    #[inline]
    fn leaf(&mut self, ids: &[u32], squared: &[f32]) {
        for (&id, &d2) in ids.iter().zip(squared) {
            if d2 <= self.bound && Some(id) != self.exclude {
                self.ids.push(id);
            }
        }
    }
}

/// Candidate search radius that covers every pair whose exact distance is
/// within `cutoff`, however the f32 image arithmetic rounds.
fn inflated(cutoff: f32) -> f32 {
    let wide = cutoff * (1.0 + 1.0e-4) + 1.0e-4;
    wide * wide
}

/// A reusable radius collector for one tree and cutoff.
pub(super) struct RadiusCollector {
    cutoff: f32,
    squared_cutoff: f32,
    images: Vec<[f64; 3]>,
    scratch: Scratch,
}

impl KdTree<'_> {
    pub(super) fn radius_collector(&self, cutoff: f32) -> Result<RadiusCollector, SpatialError> {
        let mut images = Vec::new();
        if let Some(periodic) = self.periodic {
            let limits = super::periodic_image_limits(&periodic, cutoff * (1.0 + 1.0e-4) + 1.0e-4)?;
            super::validate_image_budget(limits, self.periodic_options)?;
            self.image_translations(limits, &mut images)?;
        }
        Ok(RadiusCollector {
            cutoff,
            squared_cutoff: cutoff * cutoff,
            images,
            scratch: Scratch::default(),
        })
    }
}

impl RadiusCollector {
    /// Emits `(target, squared distance)` for every target within the cutoff
    /// of `atom`, excluding `atom` itself. Periodic results are unique per
    /// target, ordered by atom, with the exact minimum-image distance.
    pub(super) fn collect(
        &mut self,
        tree: &KdTree<'_>,
        atom: u32,
        mut emit: impl FnMut(u32, f32),
    ) -> Result<(), SpatialError> {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let Some(position) = tree.positions.get(index).copied() else {
            return Err(SpatialError::AtomOutOfBounds(atom));
        };
        if !finite(position) {
            return Ok(());
        }
        let Some(periodic) = tree.periodic else {
            let mut sink = RadiusSink {
                bound: self.squared_cutoff,
                exclude: Some(atom),
                emit: &mut emit,
            };
            tree.traverse(position, &mut self.scratch.stack, &mut sink);
            return Ok(());
        };
        let candidates = &mut self.scratch.candidates;
        candidates.clear();
        let mut sink = CandidateSink {
            bound: inflated(self.cutoff),
            exclude: Some(atom),
            ids: &mut *candidates,
        };
        tree.search_images(position, &self.images, &mut self.scratch.stack, &mut sink);
        candidates.sort_unstable();
        candidates.dedup();
        for &target in candidates.iter() {
            let Some(other) = tree.positions.get(target as usize).copied() else {
                continue;
            };
            let squared = periodic.distance_squared(position, other);
            if squared <= self.squared_cutoff {
                emit(target, squared);
            }
        }
        Ok(())
    }
}

/// Reusable nearest-neighbour searcher owning its heap and traversal scratch.
pub(super) struct KnnSearcher {
    pub(super) heap: KnnHeap,
    scratch: Scratch,
    images: Vec<[f64; 3]>,
    limit: f32,
}

impl KdTree<'_> {
    pub(super) fn knn_searcher(
        &self,
        capacity: usize,
        max_distance: Option<f32>,
    ) -> Result<KnnSearcher, SpatialError> {
        let mut images = Vec::new();
        let mut limit = f32::INFINITY;
        if let Some(periodic) = self.periodic {
            let limits = match max_distance {
                Some(distance) => {
                    super::periodic_image_limits(&periodic, distance * (1.0 + 1.0e-4) + 1.0e-4)?
                }
                None => NEAREST_IMAGE_LIMITS,
            };
            super::validate_image_budget(limits, self.periodic_options)?;
            self.image_translations(limits, &mut images)?;
        }
        if let Some(distance) = max_distance {
            // Candidate search is inflated for periodic rounding; the exact
            // filter in `finish` is authoritative.
            limit = if self.periodic.is_some() {
                inflated(distance)
            } else {
                distance * distance
            };
        }
        Ok(KnnSearcher {
            heap: KnnHeap::new(capacity),
            scratch: Scratch::default(),
            images,
            limit,
        })
    }

    /// Single-query nearest search used by `k_nearest`.
    pub(super) fn knn_into(
        &self,
        heap: &mut KnnHeap,
        point: [f32; 3],
        exclude: Option<u32>,
        limit: f32,
    ) -> Result<(), SpatialError> {
        let mut scratch = Scratch::default();
        let mut images = Vec::new();
        if self.periodic.is_some() {
            self.image_translations(NEAREST_IMAGE_LIMITS, &mut images)?;
            super::validate_image_budget(NEAREST_IMAGE_LIMITS, self.periodic_options)?;
        }
        heap.reset(limit, self.periodic.is_some());
        let mut sink = KnnSink { heap, exclude };
        self.search_images(point, &images, &mut scratch.stack, &mut sink);
        if let Some(periodic) = self.periodic {
            exact_periodic(
                self,
                &periodic,
                point,
                heap,
                &mut scratch.exact,
                f32::INFINITY,
            );
        }
        Ok(())
    }
}

/// Replaces image distances by exact minimum-image distances and re-ranks.
fn exact_periodic(
    tree: &KdTree<'_>,
    periodic: &crate::PeriodicBox,
    point: [f32; 3],
    heap: &mut KnnHeap,
    exact: &mut Vec<(f32, u32)>,
    limit: f32,
) {
    exact.clear();
    for &(_, atom) in heap.items() {
        if let Some(other) = tree.positions.get(atom as usize).copied() {
            let squared = periodic.distance_squared(point, other);
            if squared <= limit {
                exact.push((squared, atom));
            }
        }
    }
    heap.reset(f32::INFINITY, false);
    for &(squared, atom) in exact.iter() {
        heap.push(squared, atom);
    }
}

impl KnnSearcher {
    /// Computes the sorted neighbours of `atom` into `heap`.
    pub(super) fn search(
        &mut self,
        tree: &KdTree<'_>,
        atom: u32,
        exclude_self: bool,
        max_distance: Option<f32>,
    ) {
        let Some(point) = tree.positions.get(atom as usize).copied() else {
            self.heap.reset(self.limit, false);
            return;
        };
        self.heap.reset(self.limit, tree.periodic.is_some());
        if !finite(point) {
            return;
        }
        let exclude = exclude_self.then_some(atom);
        let mut sink = KnnSink {
            heap: &mut self.heap,
            exclude,
        };
        tree.search_images(point, &self.images, &mut self.scratch.stack, &mut sink);
        if let Some(periodic) = tree.periodic {
            let exact_limit = max_distance.map_or(f32::INFINITY, |distance| distance * distance);
            exact_periodic(
                tree,
                &periodic,
                point,
                &mut self.heap,
                &mut self.scratch.exact,
                exact_limit,
            );
        }
    }
}
