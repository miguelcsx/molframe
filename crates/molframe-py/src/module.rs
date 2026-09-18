use crate::api::{
    a2m_match_columns, a3m_match_columns, align_global_banded, align_global_matrix,
    align_local_matrix, align_mapping, align_region, align_semi_global_matrix,
    analyse_assign_chains, analyse_cad_contact_areas, analyse_cad_score, analyse_ce_align,
    analyse_ce_alignments, analyse_contact_map_similarity, analyse_dockq,
    analyse_equivalent_atom_mappings, analyse_gdt, analyse_gdt_ha, analyse_gdt_ts, analyse_lddt,
    analyse_ligand_symmetry_rmsd, analyse_map_chains, analyse_map_sequence_to_structure,
    analyse_qs_score, analyse_tm_score, analyse_weighted_rmsd, assign_chains, blosum62,
    cad_contact_areas, cad_score, ce_align, ce_alignments, chain_sequences, contact_map_similarity,
    decide_rmsd, dockq, equivalent_atom_mappings, gdt, gdt_ha, gdt_ts, global, global_score,
    governed_comparison_workflow, governed_interface_rmsd, governed_pocket_rmsd, interface_rmsd,
    kmer_counts, lddt, ligand_symmetry_rmsd, load_matrix, local, map_chains,
    map_sequence_to_structure, measure_mapping, minimizers, multiple_sequence_alignment,
    neighbor_joining, parse_a2m, parse_a3m, parse_clustal, parse_fasta, parse_fastq, parse_phylip,
    parse_stockholm, pocket_rmsd, progressive_msa, qs_score, read_sequence, seed_and_extend,
    semi_global, similar_kmers, tm_score, upgma, weighted_rmsd, write_a2m, write_a3m,
    write_clustal, write_fasta, write_fastq, write_phylip, write_sequence, write_stockholm,
};
use crate::atom::{PyAtom, PyAtoms};
use crate::bonds::{
    PyBondAdjacency, PyBondOrder, PyBondProvenance, PyBondRecord, PyBondTable, PyBondTableBuilder,
    PyBonds,
};
use crate::dms::{
    PyDmsBond, PyDmsCell, PyDmsFrame, PyDmsParticle, PyDmsSystem, PyDmsTopology, PyDmsVersion,
    read_dms, write_dms,
};
use crate::edit::PyCoordinateEdit;
use crate::errors::register;
use crate::geometry::{
    PyRigid, PySuperposition, angle, asphericity, asphericity_with_options, backbone_frames,
    backbone_torsions, best_fit_plane, best_fit_plane_with_options, centre_of_mass, centroid,
    circular_summary, cross, degrees, dihedral, displacement, distance, distance_matrix,
    distance_matrix_between, distance_squared, dot, gyration_axes, gyration_axes_with_options,
    helix_geometry, helix_geometry_with_options, inertia_tensor, norm, normalise, path_torsions,
    plane_deviation, plane_deviation_with_options, principal_axes, principal_axes_with_options,
    radius_of_gyration, rmsd, rmsd_flat, rmsf, rotation_mean, rotation_mean_with_options,
    superpose, superpose_with_options, symmetric, symmetric_with_options, torus_summary,
};
use crate::graph::{
    PyEdgeDirection, PyEdgeFeature, PyEdgeKind, PyGraph, PyGraphOptions, PyMissingFeaturePolicy,
    PyNodeFeature, PyNodeLevel, PySpatialBackend, build_graph,
};
use crate::hierarchy::{
    PyChain, PyChains, PyModel, PyModels, PyResidue, PyResidueAtoms, PyResidues,
};
use crate::intrinsic::{
    PyCartesianFit, PyDiffusionMap, PyPcaResult, PyPeriodicAngle, PyRotation3,
    PySurfaceGridOptions, PySurfaceMesh, PySurfaceWorkflowOptions, PySurfaceWorkflowResult,
    analyse_diffusion, analyse_pca, analyse_surface_geometry, analyse_torsion_pca, cartesian_pca,
    diffusion_map, surface_mesh, torsion_pca,
};
use crate::io::{
    PyAmbiguousResidueBoundaryPolicy, PyFormat, PyLimits, PyMissingElementPolicy, PyParseMode,
    PyPdbIdentifierNamespace, PyPdbWriteOptions, PyReadOptions, PyReadReport, PyReadScope,
};
use crate::ml::{
    PyDataset, PyDatasetEntry, PyDatasetFilter, PyDatasetSplit, PyDatasetWarning, PySplitOptions,
    PySplitRatios, PySplitStrategy, write_atom_ipc, write_atom_ipc_with_metadata,
    write_atom_parquet, write_atom_parquet_with_metadata,
};
use crate::structure::PyStructure;
use crate::trajectory::{
    PyTrajectory, PyTrajectoryFormat, PyTrajectoryUnits, PyTrajectoryWriteOptions,
};
use pyo3::prelude::*;
#[path = "core/registration/analysis.rs"]
mod analysis_registration;
#[path = "core/registration/ensemble.rs"]
mod ensemble_registration;
#[path = "core/registration/facade.rs"]
mod facade_registration;
#[path = "core/registration/fx.rs"]
mod fx_registration;
#[path = "core/registration/geometry.rs"]
mod geometry_registration;
#[path = "core/registration/governance.rs"]
mod governance_registration;
#[path = "core/registration/ic.rs"]
mod ic_registration;
#[path = "core/registration/indices.rs"]
mod index_registration;
#[path = "core/registration/io.rs"]
mod io_registration;
#[path = "core/registration/modelcif.rs"]
mod modelcif_registration;
#[path = "core/registration/namespaces.rs"]
mod namespace_registration;
#[path = "core/registration/query.rs"]
mod query_registration;
#[path = "core/registration/read.rs"]
mod read_registration;
use analysis_registration::{
    register_analysis, register_contract_classes, register_domain_classes,
};
use ensemble_registration::register_ensemble;
use facade_registration::register_facade;
use geometry_registration::register_geometry_classes;
use governance_registration::register_governance;
use ic_registration::register_ic;
use index_registration::register_indices;
use io_registration::register_io_functions;
use modelcif_registration::register_modelcif;
use query_registration::register_query;
use read_registration::read;
// No `gil_used = false` here, and the reason is a measurement rather than an
// omission.  The claim is about a free-threaded interpreter, and this repo has
// no way to produce one: the wheel is `abi3-py313`, and in the pinned pyo3
// 0.29.2 the free-threaded stable ABI starts at `abi3t-py315` — free-threaded
// CPython before 3.15 gets a version-specific build instead, so there is no
// `abi3t-py313` to ask for.  The audit the claim would rest on is done, and
// clean: nothing in the crate is `unsendable`, nothing holds a `GILProtected`,
// the two lazily created exception types are `PyOnceLock` (pyo3's
// free-threading-safe lazy global) and the single shared instance, `Columns`,
// is a frozen unit struct.  Checking with `--features pyo3/abi3t-py315` sets
// `Py_GIL_DISABLED` and compiles, but an exit-0 compile is a necessary condition
// and not a test — pyo3 does not use that cfg to bound `#[pyclass]` payloads, so
// it says nothing about the 794 hand-written kernels at runtime.  Declare it
// with a free-threaded interpreter to run against.
//
// The Rust name is `native`, not `_native`: the leading underscore belongs to the
// Python module, and putting it on the item as well would make every in-process
// caller of it trip `clippy::used_underscore_items`.  `name` keeps the exported
// module — and the `PyInit__native` symbol, and the module every `create_exception!`
// in this crate attaches to — exactly as it was.
#[pymodule]
#[pyo3(name = "_native")]
fn native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register(module)?;
    crate::adapters::register(module)?;
    register_facade(module)?;
    crate::chemistry::register(module)?;
    crate::bcif::register(module)?;
    crate::pdb_primitives::register(module)?;
    crate::audit::register(module)?;
    fx_registration::register_fx(module)?;
    crate::metadata::register(module)?;
    crate::mmtf_metadata::register(module)?;
    crate::ml_extensions::register(module)?;
    crate::crystallography::register(module)?;
    crate::spatial::register(module)?;
    register_ic(module)?;
    crate::trajectory::register(module)?;
    crate::plan::register(module)?;
    register_classes(module)?;
    register_analysis(module)?;
    register_functions(module)?;
    namespace_registration::register_late(module)
}
fn register_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::analysis::vector_field::register(module)?;
    crate::surface::components::register(module)?;
    crate::trajectory::register_interpolation(module)?;
    crate::core_diagnostic::register(module)?;
    crate::core_contract::register(module)?;
    crate::reexecution::register(module)?;
    register_contract_classes(module)?;
    crate::core_annotations::register(module)?;
    crate::core_values::register(module)?;
    crate::core_storage::register(module)?;
    crate::core::execution::register(module)?;
    crate::core::parallel::register(module)?;
    crate::core::provider::register(module)?;
    crate::core::structure_batches::register(module)?;
    crate::core_columns::register(module)?;
    crate::core_data::register(module)?;
    crate::core_encoded::register(module)?;
    crate::core_chunk_stats::register(module)?;
    crate::core_topology::register(module)?;
    crate::core_topology_root::register(module)?;
    crate::core_records::register(module)?;
    crate::core_views::register(module)?;
    module.add_class::<PyStructure>()?;
    register_indices(module)?;
    module.add_class::<PyFormat>()?;
    module.add_class::<PyParseMode>()?;
    module.add_class::<PyMissingElementPolicy>()?;
    module.add_class::<PyAmbiguousResidueBoundaryPolicy>()?;
    module.add_class::<PyLimits>()?;
    module.add_class::<PyReadScope>()?;
    module.add_class::<PyReadOptions>()?;
    module.add_class::<PyReadReport>()?;
    crate::core_io::register(module)?;
    module.add_class::<PyPdbIdentifierNamespace>()?;
    module.add_class::<PyPdbWriteOptions>()?;
    module.add("PdbOptions", module.getattr("PdbWriteOptions")?)?;
    module.add_class::<PyCoordinateEdit>()?;
    crate::core_edit::register(module)?;
    module.add_class::<PyAtoms>()?;
    module.add_class::<PyAtom>()?;
    module.add_class::<PyBonds>()?;
    module.add_class::<PyBondOrder>()?;
    module.add_class::<PyBondProvenance>()?;
    module.add_class::<PyBondRecord>()?;
    module.add_class::<PyBondTable>()?;
    module.add_class::<PyBondTableBuilder>()?;
    module.add_class::<PyBondAdjacency>()?;
    module.add_class::<PyModels>()?;
    module.add_class::<PyModel>()?;
    module.add_class::<PyChains>()?;
    module.add_class::<PyChain>()?;
    module.add_class::<PyResidues>()?;
    module.add_class::<PyResidue>()?;
    module.add_class::<PyResidueAtoms>()?;
    module.add_class::<PyPeriodicAngle>()?;
    module.add_class::<PyRotation3>()?;
    module.add_class::<PySurfaceMesh>()?;
    module.add_class::<PySurfaceWorkflowResult>()?;
    module.add_class::<PySurfaceWorkflowOptions>()?;
    module.add_class::<PySurfaceGridOptions>()?;
    crate::surface_types::register_classes(module)?;
    module.add_class::<PyCartesianFit>()?;
    module.add_class::<PyPcaResult>()?;
    module.add_class::<PyDiffusionMap>()?;
    module.add_class::<PyTrajectory>()?;
    module.add_class::<PyTrajectoryFormat>()?;
    module.add_class::<PyTrajectoryUnits>()?;
    module.add_class::<PyTrajectoryWriteOptions>()?;
    module.add_class::<PyDmsSystem>()?;
    module.add_class::<PyDmsParticle>()?;
    module.add_class::<PyDmsVersion>()?;
    module.add_class::<PyDmsBond>()?;
    module.add_class::<PyDmsCell>()?;
    module.add_class::<PyDmsFrame>()?;
    module.add_class::<PyDmsTopology>()?;
    module.add_class::<PyNodeLevel>()?;
    module.add_class::<PyEdgeKind>()?;
    module.add_class::<PyEdgeDirection>()?;
    module.add_class::<PyNodeFeature>()?;
    module.add_class::<PyEdgeFeature>()?;
    module.add_class::<PyMissingFeaturePolicy>()?;
    module.add_class::<PySpatialBackend>()?;
    module.add_class::<PyGraphOptions>()?;
    module.add_class::<PyGraph>()?;
    crate::ml::register_classes(module)?;
    module.add_class::<PyDatasetEntry>()?;
    module.add("ManifestEntry", module.getattr("DatasetEntry")?)?;
    module.add_class::<PyDatasetFilter>()?;
    module.add_class::<PySplitRatios>()?;
    module.add_class::<PySplitOptions>()?;
    module.add_class::<PySplitStrategy>()?;
    module.add_class::<PyDatasetWarning>()?;
    module.add_class::<PyDataset>()?;
    module.add_class::<PyDatasetSplit>()?;
    register_query(module)?;
    module.add_class::<PyRigid>()?;
    module.add_class::<PySuperposition>()?;
    register_geometry_classes(module)?;
    register_domain_classes(module)
}
fn register_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register_io_and_data_functions(module)?;
    register_compare_functions(module)?;
    register_sequence_functions(module)?;
    register_compare_aliases(module)?;
    register_surface_aliases(module)?;
    register_validation_aliases(module)?;
    module.setattr("AtomSelection", module.getattr("Selection")?)?;
    namespace_registration::register(module)
}

fn register_io_and_data_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(read, module)?)?;
    register_modelcif(module)?;
    register_io_functions(module)?;
    module.add_function(wrap_pyfunction!(write_atom_ipc, module)?)?;
    module.add_function(wrap_pyfunction!(write_atom_ipc_with_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(write_atom_parquet, module)?)?;
    module.add_function(wrap_pyfunction!(write_atom_parquet_with_metadata, module)?)?;
    module.add_function(wrap_pyfunction!(build_graph, module)?)?;
    module.add_function(wrap_pyfunction!(read_dms, module)?)?;
    module.add_function(wrap_pyfunction!(write_dms, module)?)?;
    module.add_function(wrap_pyfunction!(surface_mesh, module)?)?;
    crate::surface_functions::register_functions(module)?;
    module.add_function(wrap_pyfunction!(cartesian_pca, module)?)?;
    module.add_function(wrap_pyfunction!(torsion_pca, module)?)?;
    module.add_function(wrap_pyfunction!(diffusion_map, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_surface_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_pca, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_torsion_pca, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_diffusion, module)?)?;
    register_geometry_functions(module)?;
    Ok(())
}

fn register_compare_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(rmsd_flat, module)?)?;
    module.add_function(wrap_pyfunction!(superpose, module)?)?;
    module.add_function(wrap_pyfunction!(tm_score, module)?)?;
    module.add_function(wrap_pyfunction!(gdt_ts, module)?)?;
    module.add_function(wrap_pyfunction!(gdt_ha, module)?)?;
    module.add_function(wrap_pyfunction!(lddt, module)?)?;
    module.add_function(wrap_pyfunction!(gdt, module)?)?;
    module.add_function(wrap_pyfunction!(align_mapping, module)?)?;
    module.add_function(wrap_pyfunction!(measure_mapping, module)?)?;
    module.add_function(wrap_pyfunction!(decide_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(interface_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(pocket_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(ce_align, module)?)?;
    module.add_function(wrap_pyfunction!(ce_alignments, module)?)?;
    module.add_function(wrap_pyfunction!(qs_score, module)?)?;
    module.add_function(wrap_pyfunction!(dockq, module)?)?;
    module.add_function(wrap_pyfunction!(cad_score, module)?)?;
    module.add_function(wrap_pyfunction!(contact_map_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(equivalent_atom_mappings, module)?)?;
    module.add_function(wrap_pyfunction!(ligand_symmetry_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(chain_sequences, module)?)?;
    module.add_function(wrap_pyfunction!(map_chains, module)?)?;
    module.add_function(wrap_pyfunction!(assign_chains, module)?)?;
    module.add_function(wrap_pyfunction!(map_sequence_to_structure, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_lddt, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_tm_score, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_gdt_ts, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_gdt_ha, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_dockq, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_qs_score, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_gdt, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_ce_align, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_ce_alignments, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_map_chains, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_assign_chains, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_map_sequence_to_structure, module)?)?;
    module.add_function(wrap_pyfunction!(weighted_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_weighted_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(governed_comparison_workflow, module)?)?;
    module.add_function(wrap_pyfunction!(governed_interface_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(governed_pocket_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_cad_score, module)?)?;
    module.add_function(wrap_pyfunction!(cad_contact_areas, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_cad_contact_areas, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_contact_map_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_equivalent_atom_mappings, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_ligand_symmetry_rmsd, module)?)?;
    Ok(())
}

fn register_sequence_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(a3m_match_columns, module)?)?;
    module.add_function(wrap_pyfunction!(a2m_match_columns, module)?)?;
    module.add_function(wrap_pyfunction!(parse_a2m, module)?)?;
    module.add_function(wrap_pyfunction!(write_a2m, module)?)?;
    module.add_function(wrap_pyfunction!(parse_fastq, module)?)?;
    module.add_function(wrap_pyfunction!(write_fastq, module)?)?;
    module.add_function(wrap_pyfunction!(load_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(similar_kmers, module)?)?;
    module.add_function(wrap_pyfunction!(align_region, module)?)?;
    module.add_function(wrap_pyfunction!(align_global_banded, module)?)?;
    module.add_function(wrap_pyfunction!(global_score, module)?)?;
    module.add_function(wrap_pyfunction!(align_global_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(align_local_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(align_semi_global_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(blosum62, module)?)?;
    module.add_function(wrap_pyfunction!(kmer_counts, module)?)?;
    module.add_function(wrap_pyfunction!(minimizers, module)?)?;
    module.add_function(wrap_pyfunction!(multiple_sequence_alignment, module)?)?;
    module.add_function(wrap_pyfunction!(progressive_msa, module)?)?;
    module.add_function(wrap_pyfunction!(neighbor_joining, module)?)?;
    module.add_function(wrap_pyfunction!(parse_a3m, module)?)?;
    module.add_function(wrap_pyfunction!(parse_clustal, module)?)?;
    module.add_function(wrap_pyfunction!(parse_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(parse_phylip, module)?)?;
    module.add_function(wrap_pyfunction!(parse_stockholm, module)?)?;
    module.add_function(wrap_pyfunction!(seed_and_extend, module)?)?;
    module.add_function(wrap_pyfunction!(upgma, module)?)?;
    module.add_function(wrap_pyfunction!(write_clustal, module)?)?;
    module.add_function(wrap_pyfunction!(write_a3m, module)?)?;
    module.add_function(wrap_pyfunction!(write_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(write_phylip, module)?)?;
    module.add_function(wrap_pyfunction!(write_stockholm, module)?)?;
    module.add_function(wrap_pyfunction!(read_sequence, module)?)?;
    module.add_function(wrap_pyfunction!(write_sequence, module)?)?;
    module.add_function(wrap_pyfunction!(global, module)?)?;
    module.add_function(wrap_pyfunction!(local, module)?)?;
    module.add_function(wrap_pyfunction!(semi_global, module)?)?;
    Ok(())
}

fn register_surface_aliases(module: &Bound<'_, PyModule>) -> PyResult<()> {
    const ALIASES: &[(&str, &str)] = &[
        ("cavities_with_options", "cavities"),
        ("governed_surface_geometry", "analyse_surface_geometry"),
    ];
    for (alias, source) in ALIASES {
        module.add(*alias, module.getattr(*source)?)?;
    }
    Ok(())
}

fn register_validation_aliases(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("classify_ramachandran", module.getattr("classify")?)?;
    Ok(())
}

fn register_compare_aliases(module: &Bound<'_, PyModule>) -> PyResult<()> {
    const ALIASES: &[(&str, &str)] = &[
        ("gdt_with_cutoffs", "gdt"),
        ("lddt_with_options", "lddt"),
        ("qs_score_in_namespace", "qs_score"),
        ("dockq_in_namespace", "dockq"),
        ("governed_assign_chains", "analyse_assign_chains"),
        ("governed_cad_contact_areas", "analyse_cad_contact_areas"),
        ("governed_cad_score", "analyse_cad_score"),
        ("governed_ce_align", "analyse_ce_align"),
        ("governed_ce_alignments", "analyse_ce_alignments"),
        (
            "governed_contact_map_similarity",
            "analyse_contact_map_similarity",
        ),
        ("governed_dockq", "analyse_dockq"),
        (
            "governed_equivalent_atom_mappings",
            "analyse_equivalent_atom_mappings",
        ),
        ("governed_gdt_ha", "analyse_gdt_ha"),
        ("governed_gdt_ts", "analyse_gdt_ts"),
        ("governed_gdt_with_cutoffs", "analyse_gdt"),
        ("governed_lddt", "analyse_lddt"),
        (
            "governed_ligand_symmetry_rmsd",
            "analyse_ligand_symmetry_rmsd",
        ),
        ("governed_map_chains", "analyse_map_chains"),
        (
            "governed_map_sequence_to_structure",
            "analyse_map_sequence_to_structure",
        ),
        ("governed_qs_score", "analyse_qs_score"),
        ("governed_tm_score", "analyse_tm_score"),
        ("governed_weighted_rmsd", "analyse_weighted_rmsd"),
    ];
    for (alias, source) in ALIASES {
        module.add(*alias, module.getattr(*source)?)?;
    }
    Ok(())
}

fn register_geometry_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(distance_matrix, module)?)?;
    module.add_function(wrap_pyfunction!(distance_matrix_between, module)?)?;
    module.add_function(wrap_pyfunction!(rmsf, module)?)?;
    module.add_function(wrap_pyfunction!(distance, module)?)?;
    module.add_function(wrap_pyfunction!(distance_squared, module)?)?;
    module.add_function(wrap_pyfunction!(displacement, module)?)?;
    module.add_function(wrap_pyfunction!(dot, module)?)?;
    module.add_function(wrap_pyfunction!(cross, module)?)?;
    module.add_function(wrap_pyfunction!(norm, module)?)?;
    module.add_function(wrap_pyfunction!(normalise, module)?)?;
    module.add_function(wrap_pyfunction!(angle, module)?)?;
    module.add_function(wrap_pyfunction!(dihedral, module)?)?;
    module.add_function(wrap_pyfunction!(degrees, module)?)?;
    module.add_function(wrap_pyfunction!(centroid, module)?)?;
    module.add_function(wrap_pyfunction!(centre_of_mass, module)?)?;
    module.add_function(wrap_pyfunction!(radius_of_gyration, module)?)?;
    module.add_function(wrap_pyfunction!(inertia_tensor, module)?)?;
    module.add_function(wrap_pyfunction!(principal_axes, module)?)?;
    module.add_function(wrap_pyfunction!(asphericity, module)?)?;
    module.add_function(wrap_pyfunction!(asphericity_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(gyration_axes, module)?)?;
    module.add_function(wrap_pyfunction!(gyration_axes_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(best_fit_plane, module)?)?;
    module.add_function(wrap_pyfunction!(best_fit_plane_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(plane_deviation, module)?)?;
    module.add_function(wrap_pyfunction!(plane_deviation_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(principal_axes_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(circular_summary, module)?)?;
    module.add_function(wrap_pyfunction!(torus_summary, module)?)?;
    module.add_function(wrap_pyfunction!(rotation_mean, module)?)?;
    module.add_function(wrap_pyfunction!(rotation_mean_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(path_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(backbone_frames, module)?)?;
    module.add_function(wrap_pyfunction!(helix_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(helix_geometry_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(symmetric, module)?)?;
    module.add_function(wrap_pyfunction!(symmetric_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(backbone_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(superpose_with_options, module)?)?;
    Ok(())
}

#[cfg(test)]
#[path = "module_tests.rs"]
mod tests;
