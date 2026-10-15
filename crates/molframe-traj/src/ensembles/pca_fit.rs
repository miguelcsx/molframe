//! Parallel construction of the row-major PCA observation matrix.
//!
//! Frames are fitted in fixed blocks on the shared worker pool; the ordered
//! consumer appends each block to the single output buffer and accumulates the
//! column sums in frame order, so the buffer and the mean are bit-identical at
//! any worker count and equal to a sequential fit.
//!
//! Time `O(n * N)` for `n` frames of `N` atoms plus the Procrustes mean when
//! requested; memory one `n x 3N` `f64` buffer plus the in-flight blocks.

use crate::ensemble::{EnsembleGeometryError, generalized_procrustes_mean_source};
use crate::frame_view::FrameSource;
use crate::numeric::f64_from_usize;
use crate::pca::CartesianFit;
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockExecutionError, BlockPlan, try_for_each_block_in};
use molframe_geom::superpose;
use std::borrow::Cow;
use std::mem::size_of;

/// Frames fitted per block.
const FRAMES_PER_BLOCK: usize = 32;

/// Row-major observations and the per-feature mean.
pub(crate) struct Fitted {
    pub(crate) data: Vec<f64>,
    pub(crate) mean: Vec<f64>,
}

pub(crate) fn fitted_observations<S: FrameSource + Sync + ?Sized>(
    frames: &S,
    fit: CartesianFit<'_>,
    features: usize,
    context: &ExecutionContext,
) -> Result<Fitted, EnsembleGeometryError> {
    let reference = match fit {
        CartesianFit::None => None,
        CartesianFit::Reference(reference) => Some(Cow::Borrowed(reference)),
        CartesianFit::IterativeMean {
            tolerance,
            max_iterations,
        } => Some(Cow::Owned(generalized_procrustes_mean_source(
            frames,
            tolerance,
            max_iterations,
        )?)),
    };
    let reference = reference.as_deref();
    let count = frames.frame_count();
    let mut data = Vec::new();
    data.try_reserve_exact(
        count
            .checked_mul(features)
            .ok_or(EnsembleGeometryError::MemoryLimit)?,
    )
    .map_err(|_| EnsembleGeometryError::MemoryLimit)?;
    let mut mean = vec![0.0; features];
    let block_bytes = FRAMES_PER_BLOCK
        .checked_mul(features)
        .and_then(|entries| entries.checked_mul(size_of::<f64>()))
        .ok_or(EnsembleGeometryError::MemoryLimit)?;
    try_for_each_block_in(
        BlockPlan::new(count, FRAMES_PER_BLOCK),
        context,
        block_bytes,
        |_, range| fit_block(frames, range, reference, features),
        |block| {
            for row in block.chunks_exact(features) {
                for (sum, value) in mean.iter_mut().zip(row) {
                    *sum += value;
                }
            }
            data.extend_from_slice(&block);
            Ok(())
        },
    )
    .map_err(|error| match error {
        BlockExecutionError::Operation(error) => error,
        BlockExecutionError::Memory(_) => EnsembleGeometryError::MemoryLimit,
        _ => EnsembleGeometryError::InvalidParameter,
    })?;
    let inverse = f64_from_usize(count)
        .ok_or(EnsembleGeometryError::InvalidParameter)?
        .recip();
    for value in &mut mean {
        *value *= inverse;
    }
    Ok(Fitted { data, mean })
}

fn fit_block<S: FrameSource + ?Sized>(
    frames: &S,
    range: std::ops::Range<usize>,
    reference: Option<&[[f32; 3]]>,
    features: usize,
) -> Result<Vec<f64>, EnsembleGeometryError> {
    let mut block = Vec::with_capacity(range.len() * features);
    for index in range {
        let frame = frames
            .frame(index)
            .ok_or(EnsembleGeometryError::DimensionMismatch)?;
        match reference {
            Some(reference) => {
                let fit = superpose(frame, reference)
                    .map_err(|_| EnsembleGeometryError::DegenerateFit)?;
                for &point in frame {
                    push_point(&mut block, fit.transform.apply(point));
                }
            }
            None => {
                for &point in frame {
                    push_point(&mut block, point);
                }
            }
        }
    }
    Ok(block)
}

fn push_point(block: &mut Vec<f64>, point: [f32; 3]) {
    block.extend(point.map(f64::from));
}

/// Per-column mean of a row-major `rows x cols` buffer, summed in row order.
pub(crate) fn column_means(data: &[f64], cols: usize) -> Result<Vec<f64>, EnsembleGeometryError> {
    let rows = data.len() / cols.max(1);
    let mut mean = vec![0.0; cols];
    for row in data.chunks_exact(cols) {
        for (sum, value) in mean.iter_mut().zip(row) {
            *sum += value;
        }
    }
    let inverse = f64_from_usize(rows)
        .ok_or(EnsembleGeometryError::InvalidParameter)?
        .recip();
    for value in &mut mean {
        *value *= inverse;
    }
    Ok(mean)
}
