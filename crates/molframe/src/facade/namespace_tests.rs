//! Each namespace re-exports the native type rather than a facade wrapper,
//! except the curated structure handles, which are the facade's own.

#[test]
fn engine_core_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::engine::core::StructureData>();
}

#[cfg(feature = "chemistry")]
#[test]
fn chemistry_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::chemistry::Component>();
}

#[cfg(feature = "geometry")]
#[test]
fn geometry_namespace_exposes_native_kernels() {
    let _ = crate::geometry::distance;
}

#[cfg(feature = "crystal")]
#[test]
fn crystallography_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::crystal::Operator>();
}

#[cfg(feature = "query")]
#[test]
fn query_namespace_exposes_native_builder() {
    let _ = crate::query::col::all();
}

#[cfg(feature = "spatial")]
#[test]
fn spatial_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::spatial::SpatialPlan>();
}

#[cfg(feature = "ic")]
#[test]
fn internal_coordinate_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::ic::Hedron>();
}

#[cfg(feature = "mmcif")]
#[test]
fn cif_namespace_exposes_native_document() {
    let _ = std::mem::size_of::<crate::formats::cif::Document>();
}

#[cfg(feature = "pdb")]
#[test]
fn pdb_namespace_exposes_native_options() {
    let _ = std::mem::size_of::<crate::formats::pdb::PdbOptions>();
}

#[cfg(feature = "modelcif")]
#[test]
fn modelcif_namespace_exposes_native_model() {
    let _ = std::mem::size_of::<crate::formats::modelcif::ModelCif>();
}

#[cfg(feature = "bcif")]
#[test]
fn binary_cif_namespace_exposes_native_document() {
    let _ = std::mem::size_of::<crate::formats::bcif::BinaryDocument>();
}

#[cfg(feature = "interop")]
#[test]
fn machine_learning_namespace_exposes_native_dataset() {
    let _ = std::mem::size_of::<crate::interop::Dataset>();
}
