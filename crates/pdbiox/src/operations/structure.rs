//! Policy-bound structure operations for the facade plan.

use super::requests::ExecutionPlanError;
use pdbiox_core::contract::{Analysis, AnalysisPolicy};
use pdbiox_core::index::ResidueIndex;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_spatial::{PeriodicBox, SpatialBackend};
use std::fmt;
use std::sync::Arc;

/// A typed structure-bound operation.
#[derive(Clone)]
pub enum StructureRequest {
    /// Canonical nucleotide base-pair detection backed by an explicit CCD provider.
    BasePairs {
        /// Component provider retained by the request.
        provider: Arc<dyn pdbiox_chem::ComponentProvider>,
        /// Hydrogen-bond and canonical-pair controls.
        options: pdbiox_analysis::BasePairOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// CCD-annotated hydrogen-bond detection.
    HydrogenBonds {
        /// Geometric hydrogen-bond options.
        options: pdbiox_analysis::HydrogenBondOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Oppositely charged atom pairs within a cutoff.
    SaltBridges {
        /// Maximum anion-cation distance.
        maximum_distance: f32,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Aromatic ring stacking.
    PiStacking {
        /// Aromatic plane and distance policy.
        options: pdbiox_analysis::PiStackingOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Cation-pi interaction detection.
    CationPi {
        /// Cation and aromatic-plane policy.
        options: pdbiox_analysis::CationPiOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Solvent-mediated hydrogen-bond bridges.
    WaterBridges {
        /// Hydrogen-bond policy used to build the bridge graph.
        options: pdbiox_analysis::WaterBridgeOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Residue-level contact map.
    ContactMap {
        /// Atom distance cutoff.
        cutoff: f32,
        /// Minimum residue separation.
        minimum_separation: u32,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Chain-interface residue detection.
    ChainInterface {
        /// First chain identifier.
        first_chain: Box<str>,
        /// Second chain identifier.
        second_chain: Box<str>,
        /// Atom distance cutoff.
        cutoff: f32,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// DSSP-compatible secondary-structure assignment.
    SecondaryStructure {
        /// Hydrogen-bond and pattern parameters.
        options: pdbiox_analysis::DsspOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Half-sphere exposure for annotated polymer residues.
    HalfSphereExposure {
        /// Neighbour radius.
        radius: f32,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Nucleic-acid torsions from the semantic backbone roles.
    NucleicTorsions {
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Gaussian-network modes over an explicit atom selection.
    GaussianNetworkModel {
        /// Selected interaction sites in the output order.
        sites: AtomSelection,
        /// Contact and eigensolver controls.
        options: pdbiox_analysis::GnmOptions,
        /// Whether contacts use the structure unit cell.
        periodic: bool,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Steric-clash validation.
    #[cfg(feature = "validate")]
    Clashes {
        /// VDW overlap tolerance.
        tolerance: f32,
        /// Named radius set.
        radii: pdbiox_chem::RadiusSet,
        /// Spatial implementation.
        backend: SpatialBackend,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Covalent bond-length deviation validation.
    #[cfg(feature = "validate")]
    BondLengthDeviations {
        /// Accepted absolute deviation.
        tolerance: f32,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Cis-peptide validation.
    #[cfg(feature = "validate")]
    CisPeptides {
        /// Threshold around the trans configuration.
        threshold_degrees: f64,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Aromatic-ring planarity validation.
    #[cfg(feature = "validate")]
    Planarity {
        /// Plane-fit and deviation policy.
        options: pdbiox_validate::PlanarityOptions,
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Occupancy and B-factor quality flags.
    #[cfg(feature = "validate")]
    Quality {
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Over-coordination validation.
    #[cfg(feature = "validate")]
    Valence {
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
    /// Polymer sequence completeness validation.
    #[cfg(feature = "validate")]
    Completeness {
        /// Data and model policy.
        policy: AnalysisPolicy,
    },
}

impl fmt::Debug for StructureRequest {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::BasePairs { .. } => "BasePairs",
            Self::HydrogenBonds { .. } => "HydrogenBonds",
            Self::SaltBridges { .. } => "SaltBridges",
            Self::PiStacking { .. } => "PiStacking",
            Self::CationPi { .. } => "CationPi",
            Self::WaterBridges { .. } => "WaterBridges",
            Self::ContactMap { .. } => "ContactMap",
            Self::ChainInterface { .. } => "ChainInterface",
            Self::SecondaryStructure { .. } => "SecondaryStructure",
            Self::HalfSphereExposure { .. } => "HalfSphereExposure",
            Self::NucleicTorsions { .. } => "NucleicTorsions",
            Self::GaussianNetworkModel { .. } => "GaussianNetworkModel",
            #[cfg(feature = "validate")]
            Self::Clashes { .. } => "Clashes",
            #[cfg(feature = "validate")]
            Self::BondLengthDeviations { .. } => "BondLengthDeviations",
            #[cfg(feature = "validate")]
            Self::CisPeptides { .. } => "CisPeptides",
            #[cfg(feature = "validate")]
            Self::Planarity { .. } => "Planarity",
            #[cfg(feature = "validate")]
            Self::Quality { .. } => "Quality",
            #[cfg(feature = "validate")]
            Self::Valence { .. } => "Valence",
            #[cfg(feature = "validate")]
            Self::Completeness { .. } => "Completeness",
        };
        formatter.write_str(name)
    }
}

/// A typed result from a [`StructureRequest`].
#[derive(Clone, Debug)]
pub enum StructureValue {
    /// Canonical base-pair records with analysis metadata.
    BasePairs(Analysis<Vec<pdbiox_analysis::BasePair>>),
    /// Hydrogen-bond records with analysis metadata.
    HydrogenBonds(Analysis<Vec<pdbiox_analysis::HydrogenBond>>),
    /// Salt-bridge records with analysis metadata.
    SaltBridges(Analysis<Vec<pdbiox_analysis::SaltBridge>>),
    /// Aromatic stacking records with analysis metadata.
    PiStacking(Analysis<Vec<pdbiox_analysis::PiStacking>>),
    /// Cation-pi records with analysis metadata.
    CationPi(Analysis<Vec<pdbiox_analysis::CationPi>>),
    /// Water-bridge records with analysis metadata.
    WaterBridges(Analysis<Vec<pdbiox_analysis::WaterBridge>>),
    /// Residue contact map with analysis metadata.
    ContactMap(Analysis<pdbiox_analysis::ContactMap>),
    /// Interface residue indices with analysis metadata.
    ChainInterface(Analysis<Vec<ResidueIndex>>),
    /// Secondary-structure records with analysis metadata.
    SecondaryStructure(Analysis<Vec<pdbiox_analysis::SseRecord>>),
    /// Half-sphere exposure records with analysis metadata.
    HalfSphereExposure(Analysis<Vec<pdbiox_analysis::HalfSphereExposure>>),
    /// Nucleic-acid torsion records with analysis metadata.
    NucleicTorsions(Analysis<Vec<pdbiox_analysis::NucleicTorsions>>),
    /// Gaussian-network modes with analysis metadata.
    GaussianNetworkModel(Analysis<pdbiox_analysis::GaussianNetworkModel>),
    /// Steric-clash records with analysis metadata.
    #[cfg(feature = "validate")]
    Clashes(Analysis<Vec<pdbiox_validate::Clash>>),
    /// Bond-length deviation records with analysis metadata.
    #[cfg(feature = "validate")]
    BondLengthDeviations(Analysis<Vec<pdbiox_validate::BondDeviation>>),
    /// Cis-peptide records with analysis metadata.
    #[cfg(feature = "validate")]
    CisPeptides(Analysis<Vec<pdbiox_validate::CisPeptide>>),
    /// Planarity flags with analysis metadata.
    #[cfg(feature = "validate")]
    Planarity(Analysis<Vec<pdbiox_validate::PlanarityFlag>>),
    /// Quality flags with analysis metadata.
    #[cfg(feature = "validate")]
    Quality(Analysis<Vec<pdbiox_validate::QualityFlag>>),
    /// Valence errors with analysis metadata.
    #[cfg(feature = "validate")]
    Valence(Analysis<Vec<pdbiox_validate::ValenceError>>),
    /// Completeness records with analysis metadata.
    #[cfg(feature = "validate")]
    Completeness(Analysis<Vec<pdbiox_validate::ChainCompleteness>>),
}

pub(crate) fn execute(
    request: &StructureRequest,
    structure: &Structure,
) -> Result<StructureValue, ExecutionPlanError> {
    if let Some(value) = execute_analysis(request, structure)? {
        return Ok(value);
    }
    #[cfg(feature = "validate")]
    if let Some(value) = execute_validation(request, structure)? {
        return Ok(value);
    }
    Err(ExecutionPlanError::Governed(
        "unsupported structure operation".into(),
    ))
}

fn execute_analysis(
    request: &StructureRequest,
    structure: &Structure,
) -> Result<Option<StructureValue>, ExecutionPlanError> {
    if let Some(value) = execute_interactions(request, structure)? {
        return Ok(Some(value));
    }
    match request {
        StructureRequest::SecondaryStructure { options, policy } => {
            let kernel = pdbiox_analysis::secondary_structure_kernel(options);
            Ok(Some(StructureValue::SecondaryStructure(run(
                structure, policy, &kernel,
            )?)))
        }
        StructureRequest::HalfSphereExposure {
            radius,
            backend,
            policy,
        } => Ok(Some(StructureValue::HalfSphereExposure(run(
            structure,
            policy,
            &pdbiox_analysis::half_sphere_exposure_kernel(*radius, *backend),
        )?))),
        StructureRequest::NucleicTorsions { policy } => {
            Ok(Some(StructureValue::NucleicTorsions(run(
                structure,
                policy,
                &pdbiox_analysis::nucleic_torsions_kernel(),
            )?)))
        }
        StructureRequest::GaussianNetworkModel {
            sites,
            options,
            periodic,
            policy,
        } => {
            let periodic_box = periodic_box(structure, *periodic)?;
            let kernel = pdbiox_analysis::gnm_kernel(sites, *options, periodic_box.as_ref());
            Ok(Some(StructureValue::GaussianNetworkModel(run(
                structure, policy, &kernel,
            )?)))
        }
        _ => Ok(None),
    }
}

fn execute_interactions(
    request: &StructureRequest,
    structure: &Structure,
) -> Result<Option<StructureValue>, ExecutionPlanError> {
    match request {
        StructureRequest::BasePairs {
            provider,
            options,
            policy,
        } => Ok(Some(StructureValue::BasePairs(run(
            structure,
            policy,
            &pdbiox_analysis::base_pairs_kernel(provider.as_ref(), *options),
        )?))),
        StructureRequest::HydrogenBonds { options, policy } => {
            Ok(Some(StructureValue::HydrogenBonds(run(
                structure,
                policy,
                &pdbiox_analysis::hydrogen_bonds_kernel(*options),
            )?)))
        }
        StructureRequest::SaltBridges {
            maximum_distance,
            backend,
            policy,
        } => Ok(Some(StructureValue::SaltBridges(run(
            structure,
            policy,
            &pdbiox_analysis::salt_bridges_kernel(*maximum_distance, *backend),
        )?))),
        StructureRequest::PiStacking { options, policy } => {
            Ok(Some(StructureValue::PiStacking(run(
                structure,
                policy,
                &pdbiox_analysis::pi_stacking_kernel(*options),
            )?)))
        }
        StructureRequest::CationPi { options, policy } => Ok(Some(StructureValue::CationPi(run(
            structure,
            policy,
            &pdbiox_analysis::cation_pi_kernel(*options),
        )?))),
        StructureRequest::WaterBridges { options, policy } => {
            Ok(Some(StructureValue::WaterBridges(run(
                structure,
                policy,
                &pdbiox_analysis::water_bridges_kernel(*options),
            )?)))
        }
        StructureRequest::ContactMap {
            cutoff,
            minimum_separation,
            backend,
            policy,
        } => Ok(Some(StructureValue::ContactMap(run(
            structure,
            policy,
            &pdbiox_analysis::contact_map_kernel(*cutoff, *minimum_separation, *backend),
        )?))),
        StructureRequest::ChainInterface {
            first_chain,
            second_chain,
            cutoff,
            backend,
            policy,
        } => Ok(Some(StructureValue::ChainInterface(run(
            structure,
            policy,
            &pdbiox_analysis::chain_interface_kernel(first_chain, second_chain, *cutoff, *backend),
        )?))),
        _ => Ok(None),
    }
}

fn periodic_box(
    structure: &Structure,
    periodic: bool,
) -> Result<Option<PeriodicBox>, ExecutionPlanError> {
    if !periodic {
        return Ok(None);
    }
    let cell = structure
        .data()
        .cell
        .ok_or_else(|| ExecutionPlanError::Governed("periodic GNM requires a unit cell".into()))?;
    PeriodicBox::from_cell(cell)
        .map(Some)
        .map_err(|error| ExecutionPlanError::Governed(error.to_string().into()))
}

#[cfg(feature = "validate")]
fn execute_validation(
    request: &StructureRequest,
    structure: &Structure,
) -> Result<Option<StructureValue>, ExecutionPlanError> {
    match request {
        #[cfg(feature = "validate")]
        StructureRequest::Clashes {
            tolerance,
            radii,
            backend,
            policy,
        } => Ok(Some(StructureValue::Clashes(run(
            structure,
            policy,
            &pdbiox_validate::clashes_kernel(*tolerance, *radii, *backend),
        )?))),
        #[cfg(feature = "validate")]
        StructureRequest::BondLengthDeviations { tolerance, policy } => {
            Ok(Some(StructureValue::BondLengthDeviations(run(
                structure,
                policy,
                &pdbiox_validate::bond_length_deviations_kernel(*tolerance),
            )?)))
        }
        #[cfg(feature = "validate")]
        StructureRequest::CisPeptides {
            threshold_degrees,
            policy,
        } => Ok(Some(StructureValue::CisPeptides(run(
            structure,
            policy,
            &pdbiox_validate::cis_peptides_kernel(*threshold_degrees),
        )?))),
        #[cfg(feature = "validate")]
        StructureRequest::Planarity { options, policy } => {
            Ok(Some(StructureValue::Planarity(run(
                structure,
                policy,
                &pdbiox_validate::planarity_kernel(*options),
            )?)))
        }
        #[cfg(feature = "validate")]
        StructureRequest::Quality { policy } => Ok(Some(StructureValue::Quality(run(
            structure,
            policy,
            &pdbiox_validate::quality_flags_kernel(),
        )?))),
        #[cfg(feature = "validate")]
        StructureRequest::Valence { policy } => Ok(Some(StructureValue::Valence(run(
            structure,
            policy,
            &pdbiox_validate::valence_kernel(),
        )?))),
        #[cfg(feature = "validate")]
        StructureRequest::Completeness { policy } => Ok(Some(StructureValue::Completeness(run(
            structure,
            policy,
            &pdbiox_validate::completeness_kernel(),
        )?))),
        _ => Ok(None),
    }
}

fn run<K>(
    structure: &Structure,
    policy: &AnalysisPolicy,
    kernel: &K,
) -> Result<Analysis<K::Output>, ExecutionPlanError>
where
    K: pdbiox_analysis::StructureKernel,
    K::Error: fmt::Display,
{
    pdbiox_analysis::analyse_structure(structure, policy, kernel)
        .map_err(|error| ExecutionPlanError::Governed(error.to_string().into()))
}
