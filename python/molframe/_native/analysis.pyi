"""The compiled ``molframe._native.analysis`` submodule."""

from .._facade_models import (
    Contacts,
)
from .._facade_structure_io import (
    base_pairs, cation_pi, centre_of_mass_radial_distribution, chain_interface, coordination_numbers, gaussian_network_model,
    half_sphere_exposure, hydrogen_bonds, native_contact_fraction, nucleic_torsions, pi_stacking, radial_distribution,
    residue_contact_map, salt_bridges, secondary_structure, surface_contacts, water_bridges,
)
from .._io_types import (
    BasePairError, CationPiError, DensityError, DsspBinaryError, DsspError, DynamicsError,
    EnsembleGeometryError, EnsembleSimilarityError, EnsembleStatisticsError, FragmentMappingError, GnmError, GovernedAnalysisError,
    GovernedEnsembleError, GovernedNativeError, HelicalError, HseError, HydrogenBondError, KMeansError,
    NativeError, NucleicTorsionError, PhysicalKernelError, PiStackingError, PolymerError, PoreError,
    RadialError, StandaloneAnalysisError,
)
from .._trajectory_kernels import (
    dielectric_from_dipoles, generalized_procrustes_mean, kmeans, pairwise_fitted_rmsd, pairwise_torus_distance, rmsd_to_reference,
)
from .._trajectory_runtime import (
    agglomerative_clustering, block_convergence, cluster_population_similarity, dbscan_clustering, group_coordinate_variance, harmonic_ensemble_similarity,
    medoid,
)
from ..analysis import (
    AssignSecondaryStructure, BasePairs, CationPiAnalysis, ChainInterfaceAnalysis, CoordinationNumbers, Gnm,
    HalfSphereExposureAnalysis, HydrogenBonds, Leaflets, LinearDensity, NucleicTorsionAnalysis, PiStackingAnalysis,
    PoreProfile, RadialDistribution, ResidueContacts, SaltBridges, SasaStreamError, StreamlineDirection,
    StreamlineOptions, SurfaceContacts, VectorFieldError, VectorFieldGrid, WaterBridges, integrate_streamlines,
    sasa_stream,
)
from ..analysis._analysis_models import (
    CentreGroup, Contact, ContactMap, DsspOptions, DsspSegment, Leaflet,
    LeafletOptions, LinearDensityOptions, PolymerStatistics, PoreOptions, PoreProfileOptions, PoreSample,
    RadialDistributionOptions, SecondaryStructure, SseKind,
)
from ..analysis._analysis_operations import (
    BaseFrame, CationPi, CationPiOptions, ContactTable, DielectricOptions, DielectricResult,
    HelicalOptions, HelicalParameters, HydrogenBond, HydrogenBondOptions, PiStacking, PiStackingOptions,
    SaltBridge, WaterBridge, WaterBridgeOptions, WaterDynamicsOptions, WaterLag, analyse_centre_of_mass_radial_distribution,
    analyse_dielectric_from_dipoles, analyse_helical_parameters, analyse_helical_steps, analyse_water_dynamics, atom_contacts, atom_contacts_between,
    density_map, helical_parameters, helical_steps, identify_leaflets, linear_density, parse_dssp_output,
    run_dssp, water_dynamics,
)
from ..analysis._governed import (
    analyse_chain_interface, analyse_contacts, analyse_half_sphere_exposure, analyse_nucleic_torsions,
)
from ..xtal.density import (
    DensityMap,
)

__all__: list[str]
