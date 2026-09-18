//! Declarative, bounded-memory brick production for MRC/CCP4 maps.
//!
//! Catalog lookup is constant time and stores no per-brick rows. A request
//! reads and owns only one brick plus the reader's reusable encoded-row buffer,
//! so host memory is independent of the logical map and catalog sizes.

use std::fs::File;
use std::io::{Read, Seek};
use std::path::Path;

use molframe_core::provider::{ChunkId, DatasetId};

use super::{MrcBlockOptions, MrcBlockReader, MrcError};

/// Default maximum bytes in one scalar brick payload.
pub const DEFAULT_MRC_BRICK_PAYLOAD_BYTES: usize = 64 * 1024 * 1024;
/// Default maximum payload plus reusable decode workspace.
pub const DEFAULT_MRC_BRICK_WORKING_SET_BYTES: usize = 65 * 1024 * 1024;

/// Stable provider-owned identity of one scalar brick.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MapBrickId(u64);

impl MapBrickId {
    /// Creates an identity without narrowing it to a GPU-local index.
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    /// Returns the complete global identity.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// Maximum host memory retained or produced by one synchronous request.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MrcBrickBudget {
    /// Maximum scalar payload size, excluding reusable decode workspace.
    pub max_payload_bytes: usize,
    /// Maximum payload plus reusable encoded-row workspace.
    pub max_working_set_bytes: usize,
}

impl Default for MrcBrickBudget {
    fn default() -> Self {
        Self {
            max_payload_bytes: DEFAULT_MRC_BRICK_PAYLOAD_BYTES,
            max_working_set_bytes: DEFAULT_MRC_BRICK_WORKING_SET_BYTES,
        }
    }
}

/// Immutable brick layout and dirty-tracking policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MrcBrickOptions {
    /// Target interior shape; tail bricks are smaller on positive faces.
    pub interior: [u16; 3],
    /// Symmetric clamped halo width.
    pub halo: u16,
    /// Provider generation copied into every request and payload.
    pub generation: u64,
    /// Per-request host-memory ceilings.
    pub budget: MrcBrickBudget,
}

impl Default for MrcBrickOptions {
    fn default() -> Self {
        Self {
            interior: [64; 3],
            halo: 1,
            generation: 0,
            budget: MrcBrickBudget::default(),
        }
    }
}

/// Logical address of an interior brick.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MapBrickAddress {
    /// Interior origin in full-resolution canonical voxels.
    pub origin: [u64; 3],
    /// MRC block providers expose source resolution as mip zero.
    pub mip: u16,
}

/// Stored and interior dimensions for one possibly truncated tail brick.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MapBrickShape {
    /// Stored shape including the symmetric halo.
    pub stored: [u16; 3],
    /// Logical interior shape.
    pub interior: [u16; 3],
    /// Symmetric halo width.
    pub halo: u16,
    /// Number of stored scalar values, locally addressable by `u32`.
    pub voxel_count: u32,
}

/// Metadata convertible field-for-field into a renderer brick descriptor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScalarBrickMetadata {
    /// Global brick identity.
    pub id: MapBrickId,
    /// Logical address and mip level.
    pub address: MapBrickAddress,
    /// Stored and interior shape.
    pub shape: MapBrickShape,
    /// Conservative minimum from the MRC header.
    pub min: f32,
    /// Conservative maximum from the MRC header.
    pub max: f32,
    /// Monotonic provider generation.
    pub generation: u64,
}

/// Catalog row generated arithmetically without a materialized catalog.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MrcBrickDescriptor {
    /// Owning dataset identity.
    pub dataset: DatasetId,
    /// Independently requested chunk identity.
    pub chunk: ChunkId,
    /// Sparse brick metadata.
    pub metadata: ScalarBrickMetadata,
}

/// One immutable scalar payload with a single contiguous allocation.
#[derive(Debug)]
pub struct ScalarBrickPayload {
    descriptor: MrcBrickDescriptor,
    actual_range: [f32; 2],
    values: Vec<f32>,
}

impl ScalarBrickPayload {
    /// Catalog descriptor associated with the payload.
    #[must_use]
    pub const fn descriptor(&self) -> MrcBrickDescriptor {
        self.descriptor
    }

    /// Exact finite minimum and maximum measured while producing this payload.
    #[must_use]
    pub const fn actual_range(&self) -> [f32; 2] {
        self.actual_range
    }

    /// Read-only contiguous x-fastest scalars, including halos.
    #[must_use]
    pub fn values(&self) -> &[f32] {
        &self.values
    }

    /// Exact scalar payload bytes retained by this value.
    #[must_use]
    pub fn payload_bytes(&self) -> usize {
        self.values.len() * size_of::<f32>()
    }
}

/// Typed planning and request failures for scalar bricks.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum MrcBrickError {
    /// Interior dimensions must all be nonzero.
    #[error("MRC brick interior dimensions must be nonzero")]
    EmptyInterior,
    /// Stored shape or identity arithmetic exceeded its representation.
    #[error("MRC brick {field} overflow")]
    Overflow {
        /// Arithmetic field that overflowed.
        field: &'static str,
    },
    /// The requested global identity is outside this provider.
    #[error("MRC brick identity {brick} is not in this provider")]
    UnknownBrick {
        /// Unknown global identity.
        brick: u64,
    },
    /// A payload exceeds its dedicated scalar budget.
    #[error("MRC brick payload requires {required} bytes; limit is {limit}")]
    PayloadBudget {
        /// Required scalar bytes.
        required: usize,
        /// Configured scalar ceiling.
        limit: usize,
    },
    /// Header extrema do not conservatively contain decoded values.
    #[error("MRC header range does not contain a decoded brick value")]
    HeaderRangeMismatch,
    /// The decoded payload contains a non-finite scalar.
    #[error("MRC brick contains a non-finite scalar")]
    NonFiniteValue,
    /// Bounded source read failed.
    #[error(transparent)]
    Mrc(#[from] MrcError),
}

/// Pull provider whose catalog and request cost do not scale with logical data.
#[derive(Debug)]
pub struct MrcBrickProvider<R> {
    reader: MrcBlockReader<R>,
    dataset: DatasetId,
    first_chunk: ChunkId,
    first_brick: MapBrickId,
    options: MrcBrickOptions,
    extent: [u64; 3],
    grid: [u64; 3],
    brick_count: u64,
}

impl MrcBrickProvider<File> {
    /// Opens an MRC map and reads only its fixed header.
    ///
    /// # Errors
    ///
    /// Returns typed I/O, header, budget, shape, or identity errors.
    pub fn open(
        path: impl AsRef<Path>,
        dataset: DatasetId,
        first_chunk: ChunkId,
        first_brick: MapBrickId,
        options: MrcBrickOptions,
    ) -> Result<Self, MrcBrickError> {
        let reader = MrcBlockReader::open(
            path,
            MrcBlockOptions::default().with_memory_limit(options.budget.max_working_set_bytes),
        )?;
        Self::new(reader, dataset, first_chunk, first_brick, options)
    }
}

impl<R: Read + Seek> MrcBrickProvider<R> {
    /// Builds an implicit mip-zero catalog from an already opened block reader.
    ///
    /// # Errors
    ///
    /// Rejects invalid shapes, budgets, host conversions, or global ID ranges.
    pub fn new(
        reader: MrcBlockReader<R>,
        dataset: DatasetId,
        first_chunk: ChunkId,
        first_brick: MapBrickId,
        options: MrcBrickOptions,
    ) -> Result<Self, MrcBrickError> {
        if options.interior.contains(&0) {
            return Err(MrcBrickError::EmptyInterior);
        }
        let [extent_x, extent_y, extent_z] = reader
            .descriptor()
            .dimensions
            .map(|value| u64::try_from(value).map_err(|_| overflow("logical extent")));
        let extent = [extent_x?, extent_y?, extent_z?];
        let grid = [
            extent[0].div_ceil(u64::from(options.interior[0])),
            extent[1].div_ceil(u64::from(options.interior[1])),
            extent[2].div_ceil(u64::from(options.interior[2])),
        ];
        let brick_count = checked_product(grid, "brick count")?;
        checked_identity_range(first_chunk.get(), brick_count, "chunk identity")?;
        checked_identity_range(first_brick.get(), brick_count, "brick identity")?;
        validate_maximum_shape(options)?;
        Ok(Self {
            reader,
            dataset,
            first_chunk,
            first_brick,
            options,
            extent,
            grid,
            brick_count,
        })
    }

    /// Owning dataset identity.
    #[must_use]
    pub const fn dataset(&self) -> DatasetId {
        self.dataset
    }

    /// Full-resolution logical extent without voxel materialization.
    #[must_use]
    pub const fn logical_extent(&self) -> [u64; 3] {
        self.extent
    }

    /// Number of implicit mip-zero bricks.
    #[must_use]
    pub const fn brick_count(&self) -> u64 {
        self.brick_count
    }

    /// Returns the bounded block reader without retaining catalog state.
    #[must_use]
    pub fn into_reader(self) -> MrcBlockReader<R> {
        self.reader
    }

    /// Generates one catalog descriptor in constant time.
    ///
    /// # Errors
    ///
    /// Returns [`MrcBrickError::UnknownBrick`] outside the provider range.
    pub fn descriptor(&self, id: MapBrickId) -> Result<MrcBrickDescriptor, MrcBrickError> {
        let ordinal = id
            .get()
            .checked_sub(self.first_brick.get())
            .filter(|value| *value < self.brick_count)
            .ok_or(MrcBrickError::UnknownBrick { brick: id.get() })?;
        let plane = self.grid[0]
            .checked_mul(self.grid[1])
            .ok_or_else(|| overflow("brick plane"))?;
        let index = [
            ordinal % self.grid[0],
            (ordinal / self.grid[0]) % self.grid[1],
            ordinal / plane,
        ];
        let origin = [
            index[0] * u64::from(self.options.interior[0]),
            index[1] * u64::from(self.options.interior[1]),
            index[2] * u64::from(self.options.interior[2]),
        ];
        let interior = [
            tail(self.extent[0], origin[0], self.options.interior[0])?,
            tail(self.extent[1], origin[1], self.options.interior[1])?,
            tail(self.extent[2], origin[2], self.options.interior[2])?,
        ];
        let shape = shape(interior, self.options.halo)?;
        let range = self.reader.descriptor().value_range;
        Ok(MrcBrickDescriptor {
            dataset: self.dataset,
            chunk: ChunkId::new(self.first_chunk.get() + ordinal),
            metadata: ScalarBrickMetadata {
                id,
                address: MapBrickAddress { origin, mip: 0 },
                shape,
                min: range[0],
                max: range[1],
                generation: self.options.generation,
            },
        })
    }

    /// Reads exactly one requested brick and no unrelated logical payload.
    ///
    /// # Errors
    ///
    /// Returns typed identity, budget, finite-value, header-range, or I/O errors.
    pub fn read_brick(&mut self, id: MapBrickId) -> Result<ScalarBrickPayload, MrcBrickError> {
        let descriptor = self.descriptor(id)?;
        let payload_bytes = usize::try_from(descriptor.metadata.shape.voxel_count)
            .map_err(|_| overflow("payload bytes"))?
            .checked_mul(size_of::<f32>())
            .ok_or_else(|| overflow("payload bytes"))?;
        if payload_bytes > self.options.budget.max_payload_bytes {
            return Err(MrcBrickError::PayloadBudget {
                required: payload_bytes,
                limit: self.options.budget.max_payload_bytes,
            });
        }
        let [origin_x, origin_y, origin_z] = descriptor
            .metadata
            .address
            .origin
            .map(|value| usize::try_from(value).map_err(|_| overflow("host origin")));
        let origin = [origin_x?, origin_y?, origin_z?];
        let interior = descriptor.metadata.shape.interior.map(usize::from);
        let mut values = Vec::new();
        self.reader.read_block_with_halo_into(
            origin,
            interior,
            usize::from(descriptor.metadata.shape.halo),
            &mut values,
        )?;
        let actual_range = finite_range(&values)?;
        if actual_range[0] < descriptor.metadata.min || actual_range[1] > descriptor.metadata.max {
            return Err(MrcBrickError::HeaderRangeMismatch);
        }
        Ok(ScalarBrickPayload {
            descriptor,
            actual_range,
            values,
        })
    }
}

fn validate_maximum_shape(options: MrcBrickOptions) -> Result<(), MrcBrickError> {
    let shape = shape(options.interior, options.halo)?;
    let bytes = usize::try_from(shape.voxel_count)
        .map_err(|_| overflow("payload bytes"))?
        .checked_mul(size_of::<f32>())
        .ok_or_else(|| overflow("payload bytes"))?;
    if bytes > options.budget.max_payload_bytes {
        return Err(MrcBrickError::PayloadBudget {
            required: bytes,
            limit: options.budget.max_payload_bytes,
        });
    }
    Ok(())
}

fn shape(interior: [u16; 3], halo: u16) -> Result<MapBrickShape, MrcBrickError> {
    let doubled = halo.checked_mul(2).ok_or_else(|| overflow("halo"))?;
    let stored = [
        interior[0]
            .checked_add(doubled)
            .ok_or_else(|| overflow("stored shape"))?,
        interior[1]
            .checked_add(doubled)
            .ok_or_else(|| overflow("stored shape"))?,
        interior[2]
            .checked_add(doubled)
            .ok_or_else(|| overflow("stored shape"))?,
    ];
    let count = stored
        .into_iter()
        .try_fold(1_u32, |total, value| total.checked_mul(u32::from(value)));
    Ok(MapBrickShape {
        stored,
        interior,
        halo,
        voxel_count: count.ok_or_else(|| overflow("voxel count"))?,
    })
}

fn tail(extent: u64, origin: u64, target: u16) -> Result<u16, MrcBrickError> {
    u16::try_from((extent - origin).min(u64::from(target))).map_err(|_| overflow("tail shape"))
}

fn checked_product(values: [u64; 3], field: &'static str) -> Result<u64, MrcBrickError> {
    values
        .into_iter()
        .try_fold(1_u64, u64::checked_mul)
        .ok_or_else(|| overflow(field))
}

fn checked_identity_range(
    first: u64,
    count: u64,
    field: &'static str,
) -> Result<(), MrcBrickError> {
    if count > 0 && first.checked_add(count - 1).is_none() {
        return Err(overflow(field));
    }
    Ok(())
}

fn finite_range(values: &[f32]) -> Result<[f32; 2], MrcBrickError> {
    let mut range = [f32::INFINITY, f32::NEG_INFINITY];
    for value in values {
        if !value.is_finite() {
            return Err(MrcBrickError::NonFiniteValue);
        }
        range[0] = range[0].min(*value);
        range[1] = range[1].max(*value);
    }
    Ok(range)
}

const fn overflow(field: &'static str) -> MrcBrickError {
    MrcBrickError::Overflow { field }
}

#[cfg(test)]
#[path = "mrc_brick_tests.rs"]
mod tests;
