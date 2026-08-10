//! Native-module registration and facade entry points.

use crate::analysis::{
    PyAtomDepthOptions, PyBasePair, PyBasePairOptions, PyBondDeviation, PyCartesianAxis,
    PyCationPi, PyCationPiOptions, PyCavity, PyChainCompleteness, PyChiralityFlag,
    PyChiralityIssue, PyChiralityOptions, PyChiralityReport, PyCisPeptide, PyClash, PyContact,
    PyContactMap, PyDensityGrid, PyDensityGridSpec, PyDsspOptions, PyFragmentMatch,
    PyFragmentReference, PyGaussianNetworkModel, PyGnmOptions, PyHalfSphereExposure,
    PyHydrogenBond, PyHydrogenBondOptions, PyLinearDensityBin, PyMissingResidue, PyNativeContacts,
    PyNucleicTorsions, PyPiStacking, PyPiStackingOptions, PyPlanarityFlag, PyPlanarityOptions,
    PyPolymerStatistics, PyPoreOptions, PyPoreSample, PyPucker, PyQualityFlag, PyQualityIssue,
    PyRadialBin, PyRadialOptions, PyRamachandranBasin, PyRamachandranOptions, PyRamachandranRecord,
    PyRamachandranRegion, PyReferenceAssessment, PyReferenceDistribution, PyReferenceLibrary,
    PyResidueContact, PyRotamerDefinition, PyRotamerFlag, PyRotamerOptions, PyRotamerProfile,
    PyRotamerReport, PySaltBridge, PySecondaryStructure, PySseKind, PyStackingKind,
    PyStereoConfiguration, PySurfaceAreas, PySurfaceContactOptions, PyValenceError, PyWaterBridge,
    analyse_chain_interface, analyse_contacts, analyse_half_sphere_exposure,
    analyse_nucleic_torsions, assess_bond_deviation, assess_ramachandran, atom_depths,
    buried_surface, cavities, classify_ramachandran, density_map, linear_density, map_fragments,
    polymer_statistics, pore_profile, solvent_accessible_surface, sugar_pucker,
    validate_bond_lengths, validate_cis_peptides, validate_clashes, validate_completeness,
    validate_planarity, validate_quality, validate_valence,
};
use crate::atom::{PyAtom, PyAtoms};
use crate::bonds::PyBonds;
use crate::compatibility::{PyPdbParser, PyUniverse};
use crate::contract::{
    PyAnalysis, PyAssumption, PyAssumptionSource, PyCoverage, PyDiagnostic, PyImpactEstimate,
    PyProvenance, PyStatus,
};
use crate::dms::{PyDmsParticle, PyDmsSystem};
use crate::edit::PyCoordinateEdit;
use crate::errors::{read_error, register};
use crate::geometry::{
    PyAxes, PyBackboneCoordinates, PyBackboneFrame, PyBackboneTorsions, PyCircularSummary,
    PyEigenOptions, PyHelixGeometry, PyPlane, PyRigid, PySuperposition, angle, asphericity,
    backbone_frames, best_fit_plane, centre_of_mass, centroid, circular_summary, cross, degrees,
    dihedral, displacement, distance, distance_matrix, distance_matrix_between, distance_squared,
    dot, gyration_axes, helix_geometry, inertia_tensor, norm, normalise, path_torsions,
    plane_deviation, principal_axes, radius_of_gyration, rmsd, rmsf, rotation_mean, superpose,
    torus_summary,
};
use crate::graph::{
    PyEdgeDirection, PyEdgeFeature, PyEdgeKind, PyGraph, PyGraphOptions, PyMissingFeaturePolicy,
    PyNodeFeature, PyNodeLevel, PySpatialBackend,
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
    PyPdbWriteOptions, PyReadOptions, PyReadReport, PyReadScope, read_bytes, read_with_options,
    write_bcif, write_mmcif, write_pdb,
};
use crate::ml::{
    PyDataset, PyDatasetEntry, PyDatasetSplit, PyDatasetWarning, PySplitStrategy, write_atom_ipc,
    write_atom_parquet,
};
use crate::query::{
    PyAltlocPolicy, PyAnalysisPolicy, PyAssemblyChoice, PyMissingPolicy, PyModelChoice,
    PyNamespace, PyQuery, PySelection,
};
use crate::science::{
    PyAlignment, PyAlignmentMode, PyCadContact, PyCadScore, PyCeAlignment, PyCeOptions,
    PyCeSignificanceProfile, PyChainAlternative, PyChainAssignment, PyChainMapping,
    PyChainSequence, PyContactArea, PyContactSimilarity, PyDockQ, PyDockQOptions,
    PyEmptyLddtPolicy, PyEmptyQsPolicy, PyEquivalentAtomMapping, PyFastaRecord, PyFastqRecord,
    PyLadderDirection, PyLddtOptions, PyLigandRmsd, PyLocalCad, PyMatrixProfile, PyMsaOptions,
    PyQsOptions, PyRegionOptions, PyResidueMatch, PyScoring, PySequenceDocument,
    PySequenceDocumentKind, PySequenceFormat, PySimilarKmer, PySimilarKmerOptions,
    PySubstitutionMatrix, PyTree, a2m_match_columns, a3m_match_columns, align_global_banded,
    align_global_matrix, align_local_matrix, align_region, align_semi_global_matrix,
    analyse_assign_chains, analyse_cad_contact_areas, analyse_cad_score, analyse_ce_align,
    analyse_ce_alignments, analyse_contact_map_similarity, analyse_dockq,
    analyse_equivalent_atom_mappings, analyse_gdt, analyse_gdt_ha, analyse_gdt_ts, analyse_lddt,
    analyse_ligand_symmetry_rmsd, analyse_map_chains, analyse_map_sequence_to_structure,
    analyse_qs_score, analyse_tm_score, analyse_weighted_rmsd, assign_chains, blosum62,
    cad_contact_areas, cad_score, ce_align, ce_alignments, chain_sequences, contact_map_similarity,
    dockq, equivalent_atom_mappings, gdt, gdt_ha, gdt_ts, global, kmer_counts, lddt,
    ligand_symmetry_rmsd, load_matrix, local, map_chains, map_sequence_to_structure, minimizers,
    multiple_sequence_alignment, neighbor_joining, parse_a2m, parse_a3m, parse_clustal,
    parse_fasta, parse_fastq, parse_phylip, parse_stockholm, progressive_msa, qs_score,
    read_sequence, seed_and_extend, semi_global, similar_kmers, tm_score, upgma, weighted_rmsd,
    write_a2m, write_a3m, write_clustal, write_fasta, write_fastq, write_phylip, write_sequence,
    write_stockholm,
};
use crate::structure::PyStructure;
use crate::trajectory::{
    PyTrajectory, PyTrajectoryFormat, PyTrajectoryUnits, PyTrajectoryWriteOptions,
};
use pyo3::prelude::*;
use std::path::PathBuf;

#[path = "registration/ensemble.rs"]
mod ensemble_registration;
#[path = "registration/governance.rs"]
mod governance_registration;
use ensemble_registration::register_ensemble;
use governance_registration::register_governance;

#[pyfunction]
#[pyo3(signature = (path, options=None))]
fn read(py: Python<'_>, path: PathBuf, options: Option<&PyReadOptions>) -> PyResult<PyStructure> {
    let options = options.map(|value| value.0.clone());
    py.detach(move || match options {
        Some(options) => pdbiox::read_with_options(path, &options).map(|value| value.0),
        None => pdbiox::read(path),
    })
    .map(PyStructure::new)
    .map_err(|findings| read_error(py, &findings))
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register(module)?;
    crate::chemistry::register(module)?;
    crate::crystallography::register(module)?;
    crate::spatial::register(module)?;
    register_classes(module)?;
    register_analysis(module)?;
    register_functions(module)
}

fn register_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register_contract_classes(module)?;
    module.add_class::<PyStructure>()?;
    module.add_class::<PyFormat>()?;
    module.add_class::<PyParseMode>()?;
    module.add_class::<PyMissingElementPolicy>()?;
    module.add_class::<PyAmbiguousResidueBoundaryPolicy>()?;
    module.add_class::<PyLimits>()?;
    module.add_class::<PyReadScope>()?;
    module.add_class::<PyReadOptions>()?;
    module.add_class::<PyReadReport>()?;
    module.add_class::<PyPdbWriteOptions>()?;
    module.add_class::<PyCoordinateEdit>()?;
    module.add_class::<PyAtoms>()?;
    module.add_class::<PyAtom>()?;
    module.add_class::<PyBonds>()?;
    module.add_class::<PyModels>()?;
    module.add_class::<PyModel>()?;
    module.add_class::<PyChains>()?;
    module.add_class::<PyChain>()?;
    module.add_class::<PyResidues>()?;
    module.add_class::<PyResidue>()?;
    module.add_class::<PyResidueAtoms>()?;
    module.add_class::<PyPdbParser>()?;
    module.add_class::<PyUniverse>()?;
    module.add_class::<PyPeriodicAngle>()?;
    module.add_class::<PyRotation3>()?;
    module.add_class::<PySurfaceMesh>()?;
    module.add_class::<PySurfaceWorkflowResult>()?;
    module.add_class::<PySurfaceWorkflowOptions>()?;
    module.add_class::<PySurfaceGridOptions>()?;
    module.add_class::<PyCartesianFit>()?;
    module.add_class::<PyPcaResult>()?;
    module.add_class::<PyDiffusionMap>()?;
    module.add_class::<PyTrajectory>()?;
    module.add_class::<PyTrajectoryFormat>()?;
    module.add_class::<PyTrajectoryUnits>()?;
    module.add_class::<PyTrajectoryWriteOptions>()?;
    module.add_class::<PyDmsSystem>()?;
    module.add_class::<PyDmsParticle>()?;
    module.add_class::<PyNodeLevel>()?;
    module.add_class::<PyEdgeKind>()?;
    module.add_class::<PyEdgeDirection>()?;
    module.add_class::<PyNodeFeature>()?;
    module.add_class::<PyEdgeFeature>()?;
    module.add_class::<PyMissingFeaturePolicy>()?;
    module.add_class::<PySpatialBackend>()?;
    module.add_class::<PyGraphOptions>()?;
    module.add_class::<PyGraph>()?;
    module.add_class::<PyDatasetEntry>()?;
    module.add_class::<PySplitStrategy>()?;
    module.add_class::<PyDatasetWarning>()?;
    module.add_class::<PyDataset>()?;
    module.add_class::<PyDatasetSplit>()?;
    module.add_class::<PyNamespace>()?;
    module.add_class::<PyMissingPolicy>()?;
    module.add_class::<PyModelChoice>()?;
    module.add_class::<PyAltlocPolicy>()?;
    module.add_class::<PyAssemblyChoice>()?;
    module.add_class::<PyAnalysisPolicy>()?;
    module.add_class::<PySelection>()?;
    module.add_class::<PyQuery>()?;
    module.add_class::<PyRigid>()?;
    module.add_class::<PySuperposition>()?;
    register_geometry_classes(module)?;
    register_science_classes(module)
}

fn register_science_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyScoring>()?;
    module.add_class::<PyAlignment>()?;
    module.add_class::<PyEmptyLddtPolicy>()?;
    module.add_class::<PyLddtOptions>()?;
    module.add_class::<PyCeSignificanceProfile>()?;
    module.add_class::<PyCeOptions>()?;
    module.add_class::<PyCeAlignment>()?;
    module.add_class::<PyEmptyQsPolicy>()?;
    module.add_class::<PyQsOptions>()?;
    module.add_class::<PyDockQOptions>()?;
    module.add_class::<PyDockQ>()?;
    module.add_class::<PyContactArea>()?;
    module.add_class::<PyCadContact>()?;
    module.add_class::<PyLocalCad>()?;
    module.add_class::<PyCadScore>()?;
    module.add_class::<PyContactSimilarity>()?;
    module.add_class::<PyEquivalentAtomMapping>()?;
    module.add_class::<PyLigandRmsd>()?;
    module.add_class::<PyChainSequence>()?;
    module.add_class::<PyChainMapping>()?;
    module.add_class::<PyChainAlternative>()?;
    module.add_class::<PyChainAssignment>()?;
    module.add_class::<PyResidueMatch>()?;
    module.add_class::<PyFastaRecord>()?;
    module.add_class::<PyFastqRecord>()?;
    module.add_class::<PyMatrixProfile>()?;
    module.add_class::<PySimilarKmerOptions>()?;
    module.add_class::<PySimilarKmer>()?;
    module.add_class::<PyAlignmentMode>()?;
    module.add_class::<PyRegionOptions>()?;
    module.add_class::<PyMsaOptions>()?;
    module.add_class::<PySubstitutionMatrix>()?;
    module.add_class::<PyTree>()?;
    module.add_class::<PyLadderDirection>()?;
    module.add_class::<PySequenceFormat>()?;
    module.add_class::<PySequenceDocumentKind>()?;
    module.add_class::<PySequenceDocument>()?;
    Ok(())
}

fn register_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(read, module)?)?;
    module.add_function(wrap_pyfunction!(read_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(read_bytes, module)?)?;
    module.add_function(wrap_pyfunction!(write_mmcif, module)?)?;
    module.add_function(wrap_pyfunction!(write_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(write_pdb, module)?)?;
    module.add_function(wrap_pyfunction!(write_atom_ipc, module)?)?;
    module.add_function(wrap_pyfunction!(write_atom_parquet, module)?)?;
    module.add_function(wrap_pyfunction!(surface_mesh, module)?)?;
    module.add_function(wrap_pyfunction!(cartesian_pca, module)?)?;
    module.add_function(wrap_pyfunction!(torsion_pca, module)?)?;
    module.add_function(wrap_pyfunction!(diffusion_map, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_surface_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_pca, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_torsion_pca, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_diffusion, module)?)?;
    register_geometry_functions(module)?;
    module.add_function(wrap_pyfunction!(rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(superpose, module)?)?;
    module.add_function(wrap_pyfunction!(tm_score, module)?)?;
    module.add_function(wrap_pyfunction!(gdt_ts, module)?)?;
    module.add_function(wrap_pyfunction!(gdt_ha, module)?)?;
    module.add_function(wrap_pyfunction!(lddt, module)?)?;
    module.add_function(wrap_pyfunction!(gdt, module)?)?;
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
    module.add_function(wrap_pyfunction!(analyse_cad_score, module)?)?;
    module.add_function(wrap_pyfunction!(cad_contact_areas, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_cad_contact_areas, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_contact_map_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_equivalent_atom_mappings, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_ligand_symmetry_rmsd, module)?)?;
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

fn register_contract_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStatus>()?;
    module.add_class::<PyCoverage>()?;
    module.add_class::<PyImpactEstimate>()?;
    module.add_class::<PyAssumptionSource>()?;
    module.add_class::<PyAssumption>()?;
    module.add_class::<PyDiagnostic>()?;
    module.add_class::<PyProvenance>()?;
    module.add_class::<PyAnalysis>()?;
    Ok(())
}

fn register_geometry_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAxes>()?;
    module.add_class::<PyEigenOptions>()?;
    module.add_class::<PyPlane>()?;
    module.add_class::<PyCircularSummary>()?;
    module.add_class::<PyBackboneCoordinates>()?;
    module.add_class::<PyBackboneFrame>()?;
    module.add_class::<PyBackboneTorsions>()?;
    module.add_class::<PyHelixGeometry>()?;
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
    module.add_function(wrap_pyfunction!(gyration_axes, module)?)?;
    module.add_function(wrap_pyfunction!(best_fit_plane, module)?)?;
    module.add_function(wrap_pyfunction!(plane_deviation, module)?)?;
    module.add_function(wrap_pyfunction!(circular_summary, module)?)?;
    module.add_function(wrap_pyfunction!(torus_summary, module)?)?;
    module.add_function(wrap_pyfunction!(rotation_mean, module)?)?;
    module.add_function(wrap_pyfunction!(path_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(backbone_frames, module)?)?;
    module.add_function(wrap_pyfunction!(helix_geometry, module)?)?;
    Ok(())
}

fn register_analysis(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register_extended_analysis(module)?;
    register_governance(module)?;
    register_ensemble(module)?;
    module.add_class::<PyBondDeviation>()?;
    module.add_class::<PyClash>()?;
    module.add_class::<PyMissingResidue>()?;
    module.add_class::<PyChainCompleteness>()?;
    module.add_class::<PyCisPeptide>()?;
    module.add_class::<PyPlanarityFlag>()?;
    module.add_class::<PyPlanarityOptions>()?;
    module.add_class::<PyValenceError>()?;
    module.add_class::<PyReferenceAssessment>()?;
    module.add_class::<PyContact>()?;
    module.add_class::<PySseKind>()?;
    module.add_class::<PySecondaryStructure>()?;
    module.add_class::<PyDsspOptions>()?;
    module.add_class::<PyQualityIssue>()?;
    module.add_class::<PyQualityFlag>()?;
    module.add_class::<PyRamachandranRegion>()?;
    module.add_class::<PyRamachandranRecord>()?;
    module.add_class::<PyPolymerStatistics>()?;
    module.add_class::<PyPoreOptions>()?;
    module.add_class::<PyPoreSample>()?;
    module.add_class::<PyRadialOptions>()?;
    module.add_class::<PyRadialBin>()?;
    module.add_class::<PyHydrogenBondOptions>()?;
    module.add_class::<PyHydrogenBond>()?;
    module.add_class::<PySaltBridge>()?;
    module.add_class::<PyStackingKind>()?;
    module.add_class::<PyPiStacking>()?;
    module.add_class::<PyPiStackingOptions>()?;
    module.add_class::<PyCationPi>()?;
    module.add_class::<PyCationPiOptions>()?;
    module.add_class::<PyWaterBridge>()?;
    module.add_class::<PyCartesianAxis>()?;
    module.add_class::<PyLinearDensityBin>()?;
    module.add_class::<PyDensityGridSpec>()?;
    module.add_class::<PyDensityGrid>()?;
    module.add_class::<PySurfaceAreas>()?;
    module.add_class::<PyCavity>()?;
    module.add_class::<PyAtomDepthOptions>()?;
    module.add_class::<PyReferenceDistribution>()?;
    module.add_class::<PyReferenceLibrary>()?;
    module.add_class::<PyRamachandranBasin>()?;
    module.add_class::<PyRamachandranOptions>()?;
    module.add_class::<PyStereoConfiguration>()?;
    module.add_class::<PyChiralityIssue>()?;
    module.add_class::<PyChiralityOptions>()?;
    module.add_class::<PyChiralityFlag>()?;
    module.add_class::<PyChiralityReport>()?;
    module.add_class::<PyRotamerDefinition>()?;
    module.add_class::<PyRotamerProfile>()?;
    module.add_class::<PyRotamerOptions>()?;
    module.add_class::<PyRotamerFlag>()?;
    module.add_class::<PyRotamerReport>()?;
    module.add_function(wrap_pyfunction!(polymer_statistics, module)?)?;
    module.add_function(wrap_pyfunction!(pore_profile, module)?)?;
    module.add_function(wrap_pyfunction!(linear_density, module)?)?;
    module.add_function(wrap_pyfunction!(density_map, module)?)?;
    module.add_function(wrap_pyfunction!(solvent_accessible_surface, module)?)?;
    module.add_function(wrap_pyfunction!(buried_surface, module)?)?;
    module.add_function(wrap_pyfunction!(atom_depths, module)?)?;
    module.add_function(wrap_pyfunction!(cavities, module)?)?;
    module.add_function(wrap_pyfunction!(assess_bond_deviation, module)?)?;
    module.add_function(wrap_pyfunction!(assess_ramachandran, module)?)?;
    module.add_function(wrap_pyfunction!(classify_ramachandran, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_contacts, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_chain_interface, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_half_sphere_exposure, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_nucleic_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(validate_clashes, module)?)?;
    module.add_function(wrap_pyfunction!(validate_bond_lengths, module)?)?;
    module.add_function(wrap_pyfunction!(validate_cis_peptides, module)?)?;
    module.add_function(wrap_pyfunction!(validate_planarity, module)?)?;
    module.add_function(wrap_pyfunction!(validate_quality, module)?)?;
    module.add_function(wrap_pyfunction!(validate_completeness, module)?)?;
    module.add_function(wrap_pyfunction!(validate_valence, module)?)?;
    Ok(())
}

fn register_extended_analysis(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyBasePairOptions>()?;
    module.add_class::<PyBasePair>()?;
    module.add_class::<PyResidueContact>()?;
    module.add_class::<PyContactMap>()?;
    module.add_class::<PyFragmentReference>()?;
    module.add_class::<PyFragmentMatch>()?;
    module.add_class::<PyGnmOptions>()?;
    module.add_class::<PyGaussianNetworkModel>()?;
    module.add_class::<PyHalfSphereExposure>()?;
    module.add_class::<PyNativeContacts>()?;
    module.add_class::<PyNucleicTorsions>()?;
    module.add_class::<PyPucker>()?;
    module.add_class::<PySurfaceContactOptions>()?;
    module.add_function(wrap_pyfunction!(map_fragments, module)?)?;
    module.add_function(wrap_pyfunction!(sugar_pucker, module)?)?;
    Ok(())
}
