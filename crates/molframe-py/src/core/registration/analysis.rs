use super::{register_ensemble, register_governance};
use crate::analysis::{
    PyAltlocOccupancyIssue, PyAltlocOccupancyOptions, PyAltlocOccupancyRecord,
    PyAltlocOccupancyReport, PyBFactorDistribution, PyBFactorOutlier, PyBaseFrame, PyBasePair,
    PyBasePairOptions, PyBondDeviation, PyCartesianAxis, PyCationPi, PyCationPiOptions,
    PyCcdCompletenessReport, PyCentreGroup, PyChainCompleteness, PyChiralityFlag, PyChiralityIssue,
    PyChiralityOptions, PyChiralityReport, PyCisPeptide, PyClash, PyContact, PyContactMap,
    PyContactTable, PyDensityGrid, PyDensityGridSpec, PyDielectricOptions, PyDielectricResult,
    PyDsspOptions, PyDsspSegment, PyFragmentMatch, PyFragmentReference, PyGaussianNetworkModel,
    PyGnmOptions, PyHalfSphereExposure, PyHelicalOptions, PyHelicalParameters, PyHydrogenBond,
    PyHydrogenBondOptions, PyLeaflet, PyLeafletOptions, PyLigandGeometryReport, PyLinearDensityBin,
    PyLinearDensityOptions, PyMissingResidue, PyNativeContacts, PyNucleicGeometryIssue,
    PyNucleicGeometryPolicy, PyNucleicGeometryRecord, PyNucleicTorsions, PyPiStacking,
    PyPiStackingOptions, PyPlanarityFlag, PyPlanarityOptions, PyPlaneRestraint,
    PyPlaneRestraintFlag, PyPlaneRestraintReport, PyPolymerStatistics, PyPoreOptions, PyPoreSample,
    PyPucker, PyQualityFlag, PyQualityIssue, PyRadialBin, PyRadialOptions, PyRamachandranBasin,
    PyRamachandranOptions, PyRamachandranRecord, PyRamachandranRegion, PyRealSpaceCorrelation,
    PyReferenceAngleFlag, PyReferenceAssessment, PyReferenceBondFlag, PyReferenceDistribution,
    PyReferenceGeometryOptions, PyReferenceGeometryReport, PyReferenceLibrary,
    PyResidueAtomCompleteness, PyResidueContact, PyRotamerDefinition, PyRotamerFlag,
    PyRotamerOptions, PyRotamerProfile, PyRotamerReport, PySaltBridge, PySecondaryStructure,
    PySseKind, PySseRecord, PyStackingKind, PySurfaceContactOptions, PyTlsBFactorFlag,
    PyTlsBFactorReport, PyTlsGroup, PyTlsModel, PyValenceError, PyWaterBridge,
    PyWaterBridgeOptions, PyWaterDynamicsOptions, PyWaterLag, altloc_occupancy_sums,
    analyse_altloc_occupancy, analyse_b_factor_distribution, analyse_ccd_completeness,
    analyse_centre_of_mass_radial_distribution, analyse_chain_interface, analyse_contacts,
    analyse_dielectric_from_dipoles, analyse_half_sphere_exposure, analyse_helical_parameters,
    analyse_helical_steps, analyse_nucleic_torsions, analyse_plane_restraints,
    analyse_tls_b_factor_consistency, analyse_water_dynamics, assess_bond_deviation,
    assess_ramachandran, atom_contacts, atom_contacts_between, b_factor_distribution, base_pairs,
    bond_length_deviations, cation_pi, ccd_missing_atoms, centre_of_mass_radial_distribution,
    chain_interface, chirality_outliers, cis_peptides, clashes, classify, completeness,
    coordination_numbers, density_map, dielectric_from_dipoles, gaussian_network_model,
    governed_masked_real_space_correlation, governed_real_space_map_correlation,
    governed_sampled_real_space_correlation, half_sphere_exposure, helical_parameters,
    helical_steps, hydrogen_bonds, identify_leaflets, ligand_geometry, ligand_geometry_outliers,
    linear_density, linear_density_with_options, map_fragments, masked_real_space_correlation,
    native_contact_fraction, nonplanar_aromatic_rings, nucleic_acid_geometry, nucleic_torsions,
    overvalent_atoms, parse_dssp_output, pi_stacking, plane_restraint_outliers, polymer_statistics,
    pore_profile, quality_flags, radial_distribution, ramachandran, ramachandran_outliers,
    real_space_map_correlation, reference_geometry, residue_contact_map, rotamer_outliers,
    run_dssp, salt_bridges, sampled_real_space_correlation, secondary_structure, sugar_pucker,
    surface_contacts, tls_b_factor_consistency, validate_altloc_occupancy, validate_bond_lengths,
    validate_ccd_completeness, validate_cis_peptides, validate_clashes, validate_completeness,
    validate_planarity, validate_plane_restraints, validate_quality, validate_valence,
    water_bridges, water_dynamics,
};
use crate::analysis::{
    PyAnisotropicNetworkModel, PyAnmOptions, PyNormalMode, PyNormalModeSet,
    anisotropic_network_model, read_nmd, write_nmd,
};
use crate::api::{
    PyAlignment, PyAlignmentMode, PyCadContact, PyCadScore, PyCeAlignment, PyCeOptions,
    PyCeSignificanceProfile, PyChainAlternative, PyChainAssignment, PyChainMapping,
    PyChainSequence, PyColumn, PyComparisonAlignment, PyComparisonVerdict, PyContactArea,
    PyContactSimilarity, PyDistanceMeasurement, PyDockQ, PyDockQOptions, PyEmptyLddtPolicy,
    PyEmptyQsPolicy, PyEquivalentAtomMapping, PyFastaRecord, PyFastqRecord, PyInterfaceRmsd,
    PyLadderDirection, PyLddtOptions, PyLigandRmsd, PyLocalCad, PyMatrixIdentity, PyMatrixProfile,
    PyMsaOptions, PyPocketRmsd, PyPointMapping, PyPointMatch, PyQsOptions, PyRegionOptions,
    PyResidueMatch, PyScoring, PySequenceDocument, PySequenceDocumentKind, PySequenceFormat,
    PySimilarKmer, PySimilarKmerOptions, PySubstitutionMatrix, PyTree,
};
use crate::contract::{
    PyAnalysis, PyAssumption, PyAssumptionSource, PyCoverage, PyDiagnostic, PyImpactEstimate,
    PyProvenance, PyStatus,
};
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub(crate) fn register_contract_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyStatus>()?;
    module.add_class::<PyCoverage>()?;
    module.add_class::<PyImpactEstimate>()?;
    module.add_class::<PyAssumptionSource>()?;
    module.add_class::<PyAssumption>()?;
    module.add_class::<PyDiagnostic>()?;
    module.add_class::<PyProvenance>()?;
    module.add_class::<PyAnalysis>()?;
    module.getattr("ImpactEstimate")?.setattr(
        "None",
        module.getattr("ImpactEstimate")?.getattr("NoImpact")?,
    )?;
    Ok(())
}

pub(crate) fn register_domain_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::api::register_sequence_alphabet(module)?;
    crate::api::register_sequence_kmer(module)?;
    crate::api::register_sequence_errors(module)?;
    crate::api::register_tree_types(module)?;
    module.add_class::<PyScoring>()?;
    module.add_class::<PyColumn>()?;
    module.add_class::<PyAlignment>()?;
    module.add_class::<PyEmptyLddtPolicy>()?;
    module.add_class::<PyLddtOptions>()?;
    module.add_class::<PyCeSignificanceProfile>()?;
    module.add_class::<PyCeOptions>()?;
    module.add_class::<PyCeAlignment>()?;
    module.add_class::<PyPointMatch>()?;
    module.add_class::<PyPointMapping>()?;
    module.add_class::<PyComparisonAlignment>()?;
    module.add_class::<PyDistanceMeasurement>()?;
    module.add_class::<PyComparisonVerdict>()?;
    module.add_class::<PyInterfaceRmsd>()?;
    module.add_class::<PyPocketRmsd>()?;
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
    module.add_class::<PyMatrixIdentity>()?;
    module.add_class::<PySubstitutionMatrix>()?;
    module.add_class::<PyTree>()?;
    module.add_class::<PyLadderDirection>()?;
    module.add_class::<PySequenceFormat>()?;
    module.add_class::<PySequenceDocumentKind>()?;
    module.add_class::<PySequenceDocument>()?;
    Ok(())
}

pub(crate) fn register_analysis(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register_governance(module)?;
    register_ensemble(module)?;
    register_analysis_classes(module)?;
    register_extended_analysis(module)?;
    register_domain_functions(module)
}

fn register_analysis_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAltlocOccupancyOptions>()?;
    module.add_class::<PyAltlocOccupancyIssue>()?;
    module.add_class::<PyAltlocOccupancyRecord>()?;
    module.add_class::<PyAltlocOccupancyReport>()?;
    module.add_class::<PyBFactorOutlier>()?;
    module.add_class::<PyBFactorDistribution>()?;
    module.add_class::<PyTlsModel>()?;
    module.add_class::<PyTlsGroup>()?;
    module.add_class::<PyTlsBFactorFlag>()?;
    module.add_class::<PyTlsBFactorReport>()?;
    module.add_class::<PyBaseFrame>()?;
    module.add_class::<PyHelicalOptions>()?;
    module.add_class::<PyHelicalParameters>()?;
    module.add_class::<PyDielectricOptions>()?;
    module.add_class::<PyDielectricResult>()?;
    module.add_class::<PyWaterDynamicsOptions>()?;
    module.add_class::<PyWaterLag>()?;
    module.add_class::<PyCcdCompletenessReport>()?;
    module.add_class::<PyResidueAtomCompleteness>()?;
    module.add_class::<PyPlaneRestraint>()?;
    module.add_class::<PyPlaneRestraintFlag>()?;
    module.add_class::<PyPlaneRestraintReport>()?;
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
    module.add_class::<PyContactTable>()?;
    module.add_class::<PySseKind>()?;
    module.add_class::<PySecondaryStructure>()?;
    module.add_class::<PySseRecord>()?;
    module.add_class::<PyDsspOptions>()?;
    module.add_class::<PyDsspSegment>()?;
    module.add_class::<PyQualityIssue>()?;
    module.add_class::<PyQualityFlag>()?;
    module.add_class::<PyRamachandranRegion>()?;
    module.add_class::<PyRamachandranRecord>()?;
    module.add_class::<PyPolymerStatistics>()?;
    module.add_class::<PyPoreOptions>()?;
    module.add_class::<PyPoreSample>()?;
    module.add_class::<PyRadialOptions>()?;
    module.add_class::<PyRadialBin>()?;
    module.add_class::<PyCentreGroup>()?;
    module.add_class::<PyLeafletOptions>()?;
    module.add_class::<PyLeaflet>()?;
    module.add_class::<PyHydrogenBondOptions>()?;
    module.add_class::<PyHydrogenBond>()?;
    module.add_class::<PySaltBridge>()?;
    module.add_class::<PyStackingKind>()?;
    module.add_class::<PyPiStacking>()?;
    module.add_class::<PyPiStackingOptions>()?;
    module.add_class::<PyCationPi>()?;
    module.add_class::<PyCationPiOptions>()?;
    module.add_class::<PyWaterBridge>()?;
    module.add_class::<PyWaterBridgeOptions>()?;
    module.add_class::<PyCartesianAxis>()?;
    module.add_class::<PyLinearDensityBin>()?;
    module.add_class::<PyLinearDensityOptions>()?;
    module.add_class::<PyDensityGridSpec>()?;
    module.add_class::<PyDensityGrid>()?;
    module.add_class::<PyReferenceDistribution>()?;
    module.add_class::<PyReferenceLibrary>()?;
    module.add_class::<PyRamachandranBasin>()?;
    module.add_class::<PyRamachandranOptions>()?;
    module.add_class::<PyChiralityIssue>()?;
    module.add_class::<PyChiralityOptions>()?;
    module.add_class::<PyChiralityFlag>()?;
    module.add_class::<PyChiralityReport>()?;
    module.add_class::<PyRotamerDefinition>()?;
    module.add_class::<PyRotamerProfile>()?;
    module.add_class::<PyRotamerOptions>()?;
    module.add_class::<PyRotamerFlag>()?;
    module.add_class::<PyRotamerReport>()?;
    module.add_class::<PyLigandGeometryReport>()?;
    module.add_class::<PyRealSpaceCorrelation>()?;
    module.add_class::<PyReferenceGeometryOptions>()?;
    module.add_class::<PyReferenceBondFlag>()?;
    module.add_class::<PyReferenceAngleFlag>()?;
    module.add_class::<PyReferenceGeometryReport>()?;
    module.add_class::<PyNucleicGeometryPolicy>()?;
    module.add_class::<PyNucleicGeometryIssue>()?;
    module.add_class::<PyNucleicGeometryRecord>()?;
    Ok(())
}

fn register_domain_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    register_analysis_functions(module)?;
    register_analysis_validation_functions(module)
}

fn register_analysis_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(polymer_statistics, module)?)?;
    module.add_function(wrap_pyfunction!(pore_profile, module)?)?;
    module.add_function(wrap_pyfunction!(linear_density, module)?)?;
    module.add_function(wrap_pyfunction!(linear_density_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(density_map, module)?)?;
    module.add_function(wrap_pyfunction!(assess_bond_deviation, module)?)?;
    module.add_function(wrap_pyfunction!(assess_ramachandran, module)?)?;
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
    module.add_function(wrap_pyfunction!(hydrogen_bonds, module)?)?;
    module.add_function(wrap_pyfunction!(gaussian_network_model, module)?)?;
    module.add_function(wrap_pyfunction!(anisotropic_network_model, module)?)?;
    module.add_function(wrap_pyfunction!(read_nmd, module)?)?;
    module.add_function(wrap_pyfunction!(write_nmd, module)?)?;
    module.add_function(wrap_pyfunction!(nucleic_torsions, module)?)?;
    module.add_function(wrap_pyfunction!(atom_contacts, module)?)?;
    module.add_function(wrap_pyfunction!(atom_contacts_between, module)?)?;
    module.add_function(wrap_pyfunction!(base_pairs, module)?)?;
    module.add_function(wrap_pyfunction!(cation_pi, module)?)?;
    module.add_function(wrap_pyfunction!(chain_interface, module)?)?;
    module.add_function(wrap_pyfunction!(half_sphere_exposure, module)?)?;
    module.add_function(wrap_pyfunction!(native_contact_fraction, module)?)?;
    module.add_function(wrap_pyfunction!(pi_stacking, module)?)?;
    module.add_function(wrap_pyfunction!(residue_contact_map, module)?)?;
    module.add_function(wrap_pyfunction!(salt_bridges, module)?)?;
    module.add_function(wrap_pyfunction!(secondary_structure, module)?)?;
    module.add_function(wrap_pyfunction!(surface_contacts, module)?)?;
    module.add_function(wrap_pyfunction!(water_bridges, module)?)?;
    Ok(())
}

fn register_analysis_validation_functions(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(analyse_altloc_occupancy, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_b_factor_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_ccd_completeness, module)?)?;
    module.add_function(wrap_pyfunction!(
        analyse_centre_of_mass_radial_distribution,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(analyse_dielectric_from_dipoles, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_helical_parameters, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_helical_steps, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_plane_restraints, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_tls_b_factor_consistency, module)?)?;
    module.add_function(wrap_pyfunction!(analyse_water_dynamics, module)?)?;
    module.add_function(wrap_pyfunction!(b_factor_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(dielectric_from_dipoles, module)?)?;
    module.add_function(wrap_pyfunction!(helical_parameters, module)?)?;
    module.add_function(wrap_pyfunction!(helical_steps, module)?)?;
    module.add_function(wrap_pyfunction!(tls_b_factor_consistency, module)?)?;
    module.add_function(wrap_pyfunction!(validate_altloc_occupancy, module)?)?;
    module.add_function(wrap_pyfunction!(validate_ccd_completeness, module)?)?;
    module.add_function(wrap_pyfunction!(validate_plane_restraints, module)?)?;
    module.add_function(wrap_pyfunction!(altloc_occupancy_sums, module)?)?;
    module.add_function(wrap_pyfunction!(bond_length_deviations, module)?)?;
    module.add_function(wrap_pyfunction!(ccd_missing_atoms, module)?)?;
    module.add_function(wrap_pyfunction!(chirality_outliers, module)?)?;
    module.add_function(wrap_pyfunction!(cis_peptides, module)?)?;
    module.add_function(wrap_pyfunction!(clashes, module)?)?;
    module.add_function(wrap_pyfunction!(classify, module)?)?;
    module.add_function(wrap_pyfunction!(completeness, module)?)?;
    module.add_function(wrap_pyfunction!(ligand_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(ligand_geometry_outliers, module)?)?;
    module.add_function(wrap_pyfunction!(nonplanar_aromatic_rings, module)?)?;
    module.add_function(wrap_pyfunction!(nucleic_acid_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(overvalent_atoms, module)?)?;
    module.add_function(wrap_pyfunction!(plane_restraint_outliers, module)?)?;
    module.add_function(wrap_pyfunction!(quality_flags, module)?)?;
    module.add_function(wrap_pyfunction!(ramachandran, module)?)?;
    module.add_function(wrap_pyfunction!(ramachandran_outliers, module)?)?;
    module.add_function(wrap_pyfunction!(reference_geometry, module)?)?;
    module.add_function(wrap_pyfunction!(rotamer_outliers, module)?)?;
    module.add_function(wrap_pyfunction!(real_space_map_correlation, module)?)?;
    module.add_function(wrap_pyfunction!(masked_real_space_correlation, module)?)?;
    module.add_function(wrap_pyfunction!(sampled_real_space_correlation, module)?)?;
    module.add_function(wrap_pyfunction!(
        governed_real_space_map_correlation,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        governed_masked_real_space_correlation,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(
        governed_sampled_real_space_correlation,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(water_dynamics, module)?)?;
    module.add_function(wrap_pyfunction!(radial_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(
        centre_of_mass_radial_distribution,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(coordination_numbers, module)?)?;
    module.add_function(wrap_pyfunction!(identify_leaflets, module)?)?;
    module.add_function(wrap_pyfunction!(parse_dssp_output, module)?)?;
    module.add_function(wrap_pyfunction!(run_dssp, module)?)?;
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
    module.add_class::<PyAnmOptions>()?;
    module.add_class::<PyAnisotropicNetworkModel>()?;
    module.add_class::<PyNormalMode>()?;
    module.add_class::<PyNormalModeSet>()?;
    module.add_class::<PyHalfSphereExposure>()?;
    module.add_class::<PyNativeContacts>()?;
    module.add_class::<PyNucleicTorsions>()?;
    module.add_class::<PyPucker>()?;
    module.add_class::<PySurfaceContactOptions>()?;
    module.add("PoreProfileOptions", module.getattr("PoreOptions")?)?;
    module.add(
        "RadialDistributionOptions",
        module.getattr("RadialOptions")?,
    )?;
    module.add_function(wrap_pyfunction!(map_fragments, module)?)?;
    module.add_function(wrap_pyfunction!(sugar_pucker, module)?)?;
    Ok(())
}
