//! Assemblies, symmetry and crystallographic coordinate transforms.

#![forbid(unsafe_code)]

mod affine;
#[cfg(test)]
mod affine_tests;
mod assembly;
mod assembly_spatial;
mod category_transform;
mod cell;
mod crystal;
mod crystal_images;
mod expression;
mod lower;
mod map_statistics;
mod materialize;
mod mrc;
mod mtz;
mod ncs;
mod numeric;
mod reflection;
mod reflection_cif;
mod restraints;
mod space_group;
mod symmetry;
mod symmetry_inverse;
mod view;

pub use affine::AffineTransform;
pub use assembly::{ASSEMBLIES_EXTENSION, AssemblyDef, AssemblySet, Generator, Operator};
pub use assembly_spatial::AssemblyNeighbor;
pub use cell::CellTransform;
pub use crystal::{
    CrystalNeighbor, DEFAULT_CRYSTAL_IMAGE_LIMIT, crystal_neighbors,
    crystal_neighbors_with_backend, crystal_neighbors_with_limit,
};
pub use expression::{DEFAULT_INSTANCE_LIMIT, OperExpression};
pub use lower::lower_assemblies;
pub use map_statistics::{MapHistogram, MapStatistics, MapStatisticsError};
pub use materialize::INSTANCE_ID_ANNOTATION;
pub use mrc::{DensityMap, MapBoundary, MrcError};
pub use mtz::{read_mtz, write_mtz};
pub use ncs::{
    NCS_EXTENSION, NcsAtomInstance, NcsCode, NcsExt, NcsOperator, NcsSet, NcsView, lower_ncs,
};
pub use reflection::{
    ReflectionColumn, ReflectionColumnType, ReflectionDataset, ReflectionError, ReflectionTable,
    ReflectionValue,
};
pub use reflection_cif::{lower_structure_factor_cif, write_structure_factor_cif};
pub use restraints::{
    AngleRestraint, BondRestraint, ChiralRestraint, ChiralVolumeSign, MonomerLibrary,
    MonomerLibraryReadError, MonomerRestraints, PlaneAtomRestraint, PlaneRestraint, RestraintError,
    TorsionRestraint, lower_monomer_library, read_monomer_library,
};
pub use space_group::{
    SpaceGroupSetting, space_group_by_hall, space_group_setting, space_group_settings,
};
pub use symmetry::{
    Rational, SYMMETRY_EXTENSION, SymmetryExt, SymmetryOperation, SymmetrySet, lower_symmetry,
};
pub use view::{AssemblyExt, AssemblyView, AtomInstance, ChainInstance};
