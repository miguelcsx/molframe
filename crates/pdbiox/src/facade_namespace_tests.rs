//! Each namespace re-exports the native type rather than a facade wrapper.

#[test]
fn core_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::core::Structure>();
}

#[cfg(feature = "chem")]
#[test]
fn chemistry_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::chem::Component>();
}

#[cfg(feature = "geom")]
#[test]
fn geometry_namespace_exposes_native_kernels() {
    let _ = crate::geom::distance;
}

#[cfg(feature = "xtal")]
#[test]
fn crystallography_namespace_exposes_native_types() {
    let _ = std::mem::size_of::<crate::xtal::Operator>();
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
    let _ = std::mem::size_of::<crate::cif::Document>();
}

#[cfg(feature = "pdb")]
#[test]
fn pdb_namespace_exposes_native_options() {
    let _ = std::mem::size_of::<crate::pdb::PdbOptions>();
}

#[cfg(feature = "modelcif")]
#[test]
fn modelcif_namespace_exposes_native_model() {
    let _ = std::mem::size_of::<crate::modelcif::ModelCif>();
}

#[cfg(feature = "bcif")]
#[test]
fn binary_cif_namespace_exposes_native_document() {
    let _ = std::mem::size_of::<crate::bcif::BinaryDocument>();
}

#[cfg(feature = "ml")]
#[test]
fn machine_learning_namespace_exposes_native_dataset() {
    let _ = std::mem::size_of::<crate::ml::Dataset>();
}
