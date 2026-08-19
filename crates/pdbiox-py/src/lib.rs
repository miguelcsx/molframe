//! Mechanical Python bindings for the Rust facade.
//!
//! Audited pointer use is confined to NumPy/Arrow/DLPack lifetime adapters;
//! scientific kernels remain in Rust. This ABI boundary is separate from the
//! operating-system mapping boundary in `pdbiox-mmap`.

#![deny(unsafe_op_in_unsafe_fn)]

pub(crate) mod adapters;
mod analysis;
pub(crate) mod audit;
pub(crate) mod bcif;
pub(crate) mod chem;
pub(crate) mod cif;
pub(crate) mod compare;
pub(crate) mod core;
mod fx;
pub(crate) mod geom;
pub(crate) mod ic;
pub(crate) mod ml;
pub(crate) mod modelcif;
mod module;
pub(crate) mod pdb;
mod query;
pub(crate) mod seq;
pub(crate) mod spatial;
pub(crate) mod surface;
pub(crate) mod traj;
mod validate;
pub(crate) mod xtal;

pub(crate) use adapters::compatibility;
pub(crate) use bcif::codec as bcif_codec;
pub(crate) use chem::bindings as chemistry;
pub(crate) use cif::{
    document as cif_document, lexer as cif_lexer, lower as cif_lower, parser as cif_parser,
    pdbml as cif_pdbml, rows as cif_rows, small as cif_small, write as cif_write,
};
pub(crate) use compare::difference;
pub(crate) use core::{
    annotations as core_annotations, atom, bonds, chunk_stats as core_chunk_stats,
    columns as core_columns, config, contract, contract_types as core_contract, data as core_data,
    diagnostic as core_diagnostic, edit, edit_types as core_edit, encoded as core_encoded, errors,
    facade, hierarchy, index, io, io_types as core_io, metadata, mmtf_metadata,
    plan::bindings as plan, records as core_records, reexecution, storage as core_storage,
    structure, topology as core_topology, topology_root as core_topology_root,
    values as core_values, views as core_views,
};
pub(crate) use geom::{self as geometry, intrinsic_geometry as intrinsic};
pub(crate) use ic::internal_coordinates;
pub(crate) use ml::{arrow, extensions as ml_extensions, graph};
pub(crate) use modelcif::write as modelcif_write;
pub(crate) use pdb::{headers as pdb_headers, primitives as pdb_primitives};
pub(crate) use seq::science;
pub(crate) use spatial::{index as spatial_index, periodic as spatial_periodic};
pub(crate) use surface::{functions as surface_functions, types as surface_types};
pub(crate) use traj::{self as trajectory, dms, dms_models};
pub(crate) use xtal::{crystallography, maps as xtal_maps, restraints as xtal_restraints};
