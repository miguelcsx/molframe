from os import PathLike
from typing import final
from numpy import float32, float64, int64, uint32
from numpy.typing import NDArray
from numpy.ma import MaskedArray
from . import adapters, analysis, audit, bcif, cif, chem, compare, core, geom, fx, ic, ml, modelcif, pdb, query, seq, spatial, surface, traj, validate, xtal
from .audit import AlignmentPolicy, AuditPlan, AuditReport, AuditRun, ContactDefinition, DimensionSensitivity, EquivalencePolicy, HydrogenPolicy, PeriodicPolicy, PlanError, PolicyDimension, PolicyField, PolicySpace, PolicyValue, Precision, SensitiveItem, SymmetryPolicy, Tolerance
from .bcif import BcifReader, BinaryDocument
from .cif import CifWriteError, CifWriteOptions
from .core.metadata import EntityKind, EntryMetadata, PolymerKind, ReferenceAlignment, ReferenceSequence, SEQUENCE_REFERENCES_EXTENSION, SequenceMapping, SequenceReferences
from .modelcif import MODEL_CIF_EXTENSION
from ._io_types import *
from ._provider import *
from ._native import classify_ramachandran
from .compare import GdtHa, GdtTs, Lddt, TmScore
from .core import (
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION,
    COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION, HBOND_ACCEPTOR_ANNOTATION,
    HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION, PLDDT_ANNOTATION,
    POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION, STEREO_CONFIGURATION_ANNOTATION,
    Aabb, AltId, AlgorithmId, AnalysisParameters, AnnotationColumn, AtomAnnotation, AtomAnnotations, AtomChunk, AtomChunkStats,
    AtomIndex, AtomRecord, BitVec, ChainRecord, ChainTable, ChunkBuilder, ColumnKind,
    CoordinateStore, EncodedColumn, ElementMask, EntityTable, ExtensionStore, Extremes,
    BondAdjacency, BondIndex, BondOrder, BondProvenance,
    BondRecord, BondTable, BondTableBuilder, ByteSpan, ChainIndex, CoordinateBlock,
    CoordinateGeneration, DictionaryFull, DictionaryVersion, DifferenceError, EntityIndex, ExecutionContext, Fingerprint, InstanceId, Interner, ModelIndex,
    MissingResidue, ModelTable, OptionalI32, OptionalSymbol, ParameterValue, ParentMapping, Position, Presence, ProfileId, ResidueIndex, SourceRef,
    ResidueRecord, ResidueTable, Structure, StructureData, StructureView, StructureEditor,
    CoordinateEditor, CoordinateStore, StructureDifferenceOptions, SymbolId, TARGET_CHUNK_ATOMS, Topology, ValidityMask, ValueDifference,
    AtomSelection, Class, Code, ContextItem, Diagnostics, Kind, Rendered, Severity, Strictness,
    bit_width, pack, unpack_one, write_output,
    OutputOptions, DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
)
from .analysis import (
    BasePair, BasePairOptions, CartesianAxis, CationPi, CationPiOptions,
    Contact, ContactMap, ContactTable, DensityGrid, DensityGridSpec,
    DsspOptions, FragmentMatch, FragmentReference, GaussianNetworkModel,
    GnmOptions, HalfSphereExposure, HydrogenBond, HydrogenBondOptions,
    LinearDensityBin, NativeContacts, NucleicTorsions, PiStacking,
    PiStackingOptions, PolymerStatistics, PoreOptions, PoreSample, Pucker,
    RadialBin, RadialOptions, ResidueContact, SaltBridge, SecondaryStructure,
    SseKind, StackingKind, SurfaceContactOptions, WaterBridge,
    density_map, linear_density, map_fragments,
    polymer_statistics, pore_profile, sugar_pucker,
)
from .surface import (
    AtomDepthOptions, BuriedSurface, BuriedSurfaceOp, Cavity, Sasa,
    SurfaceGridOptions, atom_depths, buried_surface, cavities,
    solvent_accessible_surface,
)
from .validate import (
    BondDeviation, ChainCompleteness, ChiralityFlag, ChiralityIssue,
    ChiralityOptions, ChiralityReport, CisPeptide, Clash, PlanarityFlag,
    PlanarityOptions, QualityFlag, QualityIssue, RamachandranBasin,
    RamachandranOptions, RamachandranRecord, RamachandranRegion,
    ReferenceAssessment, ReferenceDistribution, ReferenceLibrary,
    RotamerDefinition, RotamerFlag, RotamerOptions, RotamerProfile,
    RotamerReport, ValenceError, analyse_altloc_occupancy,
    analyse_b_factor_distribution, analyse_ccd_completeness,
    analyse_plane_restraints, analyse_tls_b_factor_consistency,
    assess_bond_deviation, assess_ramachandran, validate_altloc_occupancy,
    validate_ccd_completeness, validate_plane_restraints,
)
from .chem import AutomorphismLimit, ChemistryReport, CifProvider, Component, ComponentAtom, ComponentBond, ComponentCoverage, ComponentDictionary, ComponentKind, ComponentProvider, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE, Element, ElementProperties, EquivalenceCache, EquivalenceClasses, IonicRadius, IonicSpin, MemoryProvider, PeoeAtom, PeoeAtomType, PeoeBond, PeoeError, PeoeOptions, PeoeParameterProfile, PolymerAtomRole, PolymerLinkPolicy, PolymerLinkRule, PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, RadiusSet, RadiiSet, RadiusTable, SideChainDefinition, SideChainRoles, SmartsDataError, SmartsError, SmartsMatch, SmartsPattern, StereoConfiguration, MolAtom, MolBond, Molecule, MolVersion, MolAtomMetadata, MolBondMetadata, SdfProperty, MolRecord, Mol2AtomMetadata, Mol2BondMetadata, Mol2Section, Mol2Record, MolError, Mol2Error, apply_component_chemistry, apply_polymer_role_profile, automorphisms, component_coverage, component_peoe_charges, element_properties, equivalence_classes, ionic_radii, parse_mol_record, parse_sdf_records, parse_smarts, read_ccd, side_chain_definition, vdw_radius, write_mol, write_sdf, parse_mol2_record, write_mol2, peoe_charges
from .xtal import DEFAULT_CRYSTAL_IMAGE_LIMIT, DEFAULT_INSTANCE_LIMIT, SpaceGroup, SymmetryOperation, UnitCell, space_group_by_hall, space_group_by_number, space_group_by_symbol, space_group_setting, space_group_settings
from .xtal import ASSEMBLIES_EXTENSION, AffineTransform, AssemblyDef, AssemblyNeighbor, AssemblySet, AssemblyView, AtomInstance, CellTransform, ChainInstance, CrystalNeighbor, Generator, INSTANCE_ID_ANNOTATION, NCS_EXTENSION, NcsCode, NcsOperator, NcsSet, NcsView, OperExpression, Operator, Rational, SpaceGroupSetting, SymmetrySet, SYMMETRY_EXTENSION, collect_crystal_neighbors, lower_assemblies, lower_ncs, lower_symmetry
from .ml.graph import EdgeDirection, EdgeFeature, EdgeKind, Graph, GraphOptions, MissingFeaturePolicy, NodeFeature, NodeLevel, SpatialBackend
from .geom import Asphericity, Axes, BackboneCoordinates, BackboneFrame, BackboneResidue, BackboneTorsions, Centroid, CentreOfMass, CircularSummary, Decomposition, DistanceMatrix, DistanceMatrixBetween, DistanceMatrixOp, EigenError, EigenOptions, FluctuationError, GyrationAxes, HelixGeometry, InertiaTensor, MatrixError, PeriodicError, Plane, PrincipalAxes, RadiusOfGyration, Rigid, Rmsf, RotationError, RotationMeanOptions, RotationOptions, Superposition, SuperposeError, SuperposeOptions, TorusMetric, angle, asphericity, asphericity_with_options, backbone_frames, backbone_torsions, best_fit_plane, best_fit_plane_with_options, centre_of_mass, centroid, circular_summary, cross, degrees, dihedral, displacement, distance, distance_matrix, distance_matrix_between, distance_squared, dot, gyration_axes, gyration_axes_with_options, helix_geometry, helix_geometry_with_options, inertia_tensor, norm, normalise, path_torsions, plane_deviation, plane_deviation_with_options, principal_axes, principal_axes_with_options, radius_of_gyration, rmsd, rmsd_flat, rmsf, rotation_mean, rotation_mean_with_options, superpose, superpose_with_options, symmetric, symmetric_with_options, torus_summary
from .ml import ArrowStream, AtomArrowTable, AtomTable, BondArrowTable, BondTable, ChainArrowTable, ChainTable, DLDataType, DLDevice, DLManagedTensor, DLTensor, Dataset, DatasetEntry, DatasetError, DatasetFilter, DatasetSplit, DatasetWarning, DlpackError, DlpackTensor, ExportCost, GraphError, LoadError, ManifestEntry, MolframeExtension, ResidueArrowTable, ResidueTable, SplitOptions, SplitRatios, SplitStrategy, TableFileError, build_graph, extension_name, write_atom_ipc, write_atom_ipc_with_metadata, write_atom_parquet, write_atom_parquet_with_metadata
def graph(structure: Structure, options: object) -> object: ...
from .query import AltlocPolicy, AnalysisPolicy, AssemblyChoice, Evaluation, LogicalPlan, MissingPolicy, ModelChoice, Namespace, PhysicalQuery, SelectQuery, Query, QueryBuilder, Selection, col
from .seq import *
from .compare import *
from ._sequence_compare_types import *
from ._sequence_compare_operations import *
from .spatial import AtomsWithin, AutoBackendProfile, CellGridOptions, KdPeriodicOptions, NeighborList, NeighborListOptions, NeighborPair, NeighborPairs, NeighborSkinProfile, NeighborTable, PeriodicBox, PeriodicImage, SpatialError, SpatialOption, SpatialPlan, SpatialSearchOptions, atoms_within, atoms_within_with_options, nearest_neighbors, neighbor_pairs, neighbor_pairs_with_options
from .core.contract import Analysis, Assumption, AssumptionSource, Coverage, Diagnostic, ImpactEstimate, Provenance, Status
from .analysis import analyse_chain_interface, analyse_contacts, analyse_half_sphere_exposure, analyse_nucleic_torsions
from .surface import SurfaceWorkflowOptions, SurfaceWorkflowResult, analyse_surface_geometry
from .traj import contact_counts_stream
from .traj import CartesianFit, analyse_diffusion, analyse_pca, analyse_torsion_pca
from .validate import validate_bond_lengths, validate_cis_peptides, validate_clashes, validate_completeness, validate_planarity, validate_quality, validate_valence
from ._trajectory import DiffusionMap, DmsBond, DmsCell, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion, PcaResult, PeriodicAngle, Rotation3, SurfaceMesh, Trajectory, TrajectoryFormat, TrajectoryUnits, TrajectoryWriteOptions, write_dms
from .modelcif import GlobalMetric, LocalMetric, MetricDefinition, ModelCategory, ModelCif, ModelDescription, ModelRow, PairwiseMetric, ProtocolStep, QualityMetrics, SoftwareGroup, Target, Template, lower_model_cif
from .pdb import MMTF_METADATA_EXTENSION, MmtfEntityMetadata, MmtfGroupMetadata, MmtfMetadata, MmtfOptionalField, PDB_HEADERS_EXTENSION, PdbHeaderRecord, PdbHeaders, PdbIdentifierNamespace, PdbOptions, PdbWriteOptions, mmtf_metadata, with_mmtf_metadata, write_mmtf_with_metadata

from ._facade_models import *
from ._facade_structure_io import *
from .analysis import StreamlineDirection, StreamlineOptions, VectorFieldError, VectorFieldGrid, integrate_streamlines
from .surface import SurfaceComponent, SurfaceComponentError, SurfaceComponentFilter, filter_surface_components, surface_components
from .traj import StreamFrame, run_analysis_stream, rmsf_stream, rmsd_stream, TrajectoryInterpolation, TrajectoryInterpolationError, interpolate_trajectory_frames
