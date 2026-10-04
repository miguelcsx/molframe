//! Assemblies, symmetry and crystallographic coordinate transforms.

#![forbid(unsafe_code)]

mod affine;
#[cfg(test)]
mod affine_tests;
mod assembly;
mod assembly_spatial;
mod category_transform;
mod cell;
mod cell_reduction;
mod crystal;
mod crystal_batch;
mod crystal_images;
mod crystal_system;
mod expression;
mod grid;
mod lower;
mod map_statistics;
mod materialize;
mod mrc;
mod mtz;
mod ncs;
mod numeric;
mod reflection;
mod reflection_binning;
mod reflection_cif;
mod reflection_symmetry;
mod restraints;
mod scattering;
mod scattering_table;
mod space_group;
mod structure_factor;
mod structure_factor_model;
mod symmetry;
mod symmetry_inverse;
mod view;

pub use affine::AffineTransform;
pub use assembly::{ASSEMBLIES_EXTENSION, AssemblyDef, AssemblySet, Generator, Operator};
pub use assembly_spatial::{AssemblyBond, AssemblyNeighbor};
pub use cell::CellTransform;
pub use cell_reduction::{CellReductionError, ReducedCell, niggli_reduce};
pub use crystal::{
    CrystalNeighbor, CrystalNeighborBatch, CrystalNeighborOptions, DEFAULT_CRYSTAL_IMAGE_LIMIT,
    collect_crystal_neighbors, crystal_neighbor_batches, visit_crystal_neighbors,
};
pub use crystal_batch::{
    CrystalImageBatch, CrystalImageBatchOptions, crystal_image_batches, visit_crystal_images,
};
pub use crystal_images::CrystalImage;
pub use crystal_system::{chain_atom_ranges, crystal_contact_view, operation_motion};
pub use expression::{DEFAULT_INSTANCE_LIMIT, OperExpression};
pub use grid::{CubeAtom, CubeGrid, GridError, read_cube, read_dx};
pub use lower::lower_assemblies;
pub use map_statistics::{MapHistogram, MapStatistics, MapStatisticsError};
pub use mrc::{
    DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES, DEFAULT_MRC_BRICK_PAYLOAD_BYTES,
    DEFAULT_MRC_BRICK_WORKING_SET_BYTES, DensityMap, DensitySampler, MapBoundary, MapBrickAddress,
    MapBrickId, MapBrickShape, MrcBlockOptions, MrcBlockReader, MrcBrickBudget, MrcBrickDescriptor,
    MrcBrickError, MrcBrickOptions, MrcBrickProvider, MrcError, MrcMapDescriptor,
    ScalarBrickMetadata, ScalarBrickPayload,
};
pub use mtz::{read_mtz, write_mtz};
pub use ncs::{
    NCS_EXTENSION, NcsAtomInstance, NcsCode, NcsExt, NcsOperator, NcsSet, NcsView, lower_ncs,
};
pub use reflection::{
    ReflectionColumn, ReflectionColumnType, ReflectionDataset, ReflectionError, ReflectionTable,
    ReflectionValue,
};
pub use reflection_binning::{
    BinMethod, ReflectionBinningError, ResolutionBinner, amplitude_normalizers,
};
pub use reflection_cif::{lower_structure_factor_cif, write_structure_factor_cif};
pub use reflection_symmetry::ReflectionSymmetry;
pub use restraints::{
    AngleRestraint, BondRestraint, ChiralRestraint, ChiralVolumeSign, MonomerLibrary,
    MonomerLibraryReadError, MonomerRestraints, PlaneAtomRestraint, PlaneRestraint, RestraintError,
    TorsionRestraint, lower_monomer_library, read_monomer_library,
};
pub use scattering::GaussianFormFactor;
pub use space_group::{
    SpaceGroupSetting, space_group_by_hall, space_group_by_hermann_mauguin, space_group_setting,
    space_group_settings,
};
pub use structure_factor::{
    Complex64, Displacement, ScatteringSite, StructureFactorCalculator, StructureFactorError,
};
pub use structure_factor_model::structure_factors;
pub use symmetry::{
    Rational, SYMMETRY_EXTENSION, SymmetryExt, SymmetryOperation, SymmetrySet, lower_symmetry,
};
pub use view::{AssemblyExt, AssemblyView, AtomInstance, ChainInstance};
