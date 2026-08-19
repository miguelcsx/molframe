"""Facade model, hierarchy, chemistry, CIF, and internal-coordinate types."""

from os import PathLike
from typing import TypeAlias, final
from numpy import float32, float64, int64, uint32
from numpy.typing import NDArray
from numpy.ma import MaskedArray
from . import adapters, analysis, audit, bcif, cif, chem, compare, core, geom, fx, ic, ml, modelcif, pdb, query, seq, spatial, surface, traj, validate, xtal
from .audit import AlignmentPolicy, AuditPlan, AuditReport, AuditRun, ContactDefinition, DimensionSensitivity, EquivalencePolicy, HydrogenPolicy, PeriodicPolicy, PlanError, PolicyDimension, PolicyField, PolicySpace, PolicyValue, Precision, SensitiveItem, SymmetryPolicy, Tolerance
from .bcif import BcifReader, BinaryDocument
from .cif import CifWriteError, CifWriteOptions
from .core.metadata import EntityKind, EntryMetadata, PolymerKind, ReferenceAlignment, ReferenceSequence, SEQUENCE_REFERENCES_EXTENSION, SequenceMapping, SequenceReferences
from .modelcif import MODEL_CIF_EXTENSION
_Coordinate: TypeAlias = tuple[float, float, float]
from ._io_types import *
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
    CoordinateGeneration, DictionaryFull, DictionaryVersion, DifferenceError, EntityIndex, Fingerprint, InstanceId, Interner, ModelIndex,
    MissingResidue, ModelTable, OptionalI32, OptionalSymbol, ParameterValue, ParentMapping, Position, Presence, ProfileId, ResidueIndex, SourceRef,
    ResidueRecord, ResidueTable, Structure, StructureData, StructureView, StructureEditor,
    CoordinateEditor, CoordinateStore, StructureDifferenceOptions, SymbolId, TARGET_CHUNK_ATOMS, Topology, ValidityMask, ValueDifference,
    AtomSelection, Class, Code, ContextItem, Diagnostics, Kind, Rendered, Severity, Strictness,
    bit_width, pack, unpack_one, write_output,
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
from .analysis import AssignSecondaryStructure, BasePairs, CationPiAnalysis, ChainInterfaceAnalysis, CoordinationNumbers, Gnm, HalfSphereExposureAnalysis, HydrogenBonds, Leaflets, LinearDensity, NucleicTorsionAnalysis, PiStackingAnalysis, PoreProfile, RadialDistribution, ResidueContacts, SaltBridges, SurfaceContacts, WaterBridges
from .surface import (
    AtomDepthOptions, BuriedSurface, BuriedSurfaceOp, Cavity, Sasa,
    SurfaceGridOptions, atom_depths, buried_surface, cavities,
    solvent_accessible_surface,
)
from .validate import (
    BondDeviation, BondLengthDeviations, ChainCompleteness, ChiralityFlag,
    ChiralityIssue, ChiralityOptions, ChiralityReport, CisPeptide,
    CisPeptides, Clash, Completeness, PlanarityCheck, PlanarityFlag,
    PlanarityOptions, QualityFlag, QualityFlags, QualityIssue,
    RamachandranBasin, RamachandranOptions, RamachandranRecord,
    RamachandranRegion, ReferenceAssessment, ReferenceDistribution,
    ReferenceLibrary, RotamerDefinition, RotamerFlag, RotamerOptions,
    RotamerProfile, RotamerReport, StericClashes, Valence, ValenceError,
    analyse_altloc_occupancy, analyse_b_factor_distribution,
    analyse_ccd_completeness, analyse_plane_restraints,
    analyse_tls_b_factor_consistency, assess_bond_deviation,
    assess_ramachandran, validate_altloc_occupancy,
    validate_ccd_completeness, validate_plane_restraints,
)
from .chem import AutomorphismLimit, ChemistryReport, CifProvider, Component, ComponentAtom, ComponentBond, ComponentCoverage, ComponentDictionary, ComponentKind, ComponentProvider, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE, Element, ElementProperties, EquivalenceCache, EquivalenceClasses, InferBonds, IonicRadius, IonicSpin, MemoryProvider, PeoeAtom, PeoeAtomType, PeoeBond, PeoeError, PeoeOptions, PeoeParameterProfile, PolymerAtomRole, PolymerLinkPolicy, PolymerLinkRule, PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, RadiusSet, RadiiSet, RadiusTable, SideChainDefinition, SideChainRoles, SmartsDataError, SmartsError, SmartsMatch, SmartsPattern, StereoConfiguration, MolAtom, MolBond, Molecule, MolVersion, MolAtomMetadata, MolBondMetadata, SdfProperty, MolRecord, Mol2AtomMetadata, Mol2BondMetadata, Mol2Section, Mol2Record, MolError, Mol2Error, apply_component_chemistry, apply_polymer_role_profile, automorphisms, component_coverage, component_peoe_charges, element_properties, equivalence_classes, ionic_radii, parse_mol_record, parse_sdf_records, parse_smarts, read_ccd, side_chain_definition, vdw_radius, write_mol, write_sdf, parse_mol2_record, write_mol2, peoe_charges
from .xtal import DEFAULT_CRYSTAL_IMAGE_LIMIT, DEFAULT_INSTANCE_LIMIT, SpaceGroup, SymmetryOperation, UnitCell, space_group_by_hall, space_group_by_number, space_group_by_symbol, space_group_setting, space_group_settings
from .xtal import ASSEMBLIES_EXTENSION, AffineTransform, AssemblyDef, AssemblyNeighbor, AssemblySet, AssemblyView, AtomInstance, CellTransform, ChainInstance, CrystalNeighbor, Generator, INSTANCE_ID_ANNOTATION, NCS_EXTENSION, NcsCode, NcsOperator, NcsSet, NcsView, OperExpression, Operator, Rational, SpaceGroupSetting, SymmetrySet, SYMMETRY_EXTENSION, crystal_neighbors, crystal_neighbors_with_backend, crystal_neighbors_with_limit, lower_assemblies, lower_ncs, lower_symmetry
from .ml.graph import EdgeDirection, EdgeFeature, EdgeKind, Graph, GraphOptions, MissingFeaturePolicy, NodeFeature, NodeLevel, SpatialBackend
from .geom import Axes, BackboneCoordinates, BackboneFrame, BackboneResidue, BackboneTorsions, CircularSummary, Decomposition, DistanceMatrix, EigenError, EigenOptions, FluctuationError, HelixGeometry, MatrixError, PeriodicError, Plane, Rigid, RotationError, RotationMeanOptions, RotationOptions, Superposition, SuperposeError, SuperposeOptions, TorusMetric, angle, asphericity, asphericity_with_options, backbone_frames, backbone_torsions, best_fit_plane, best_fit_plane_with_options, centre_of_mass, centroid, circular_summary, cross, degrees, dihedral, displacement, distance, distance_matrix, distance_matrix_between, distance_squared, dot, gyration_axes, gyration_axes_with_options, helix_geometry, helix_geometry_with_options, inertia_tensor, norm, normalise, path_torsions, plane_deviation, plane_deviation_with_options, principal_axes, principal_axes_with_options, radius_of_gyration, rmsd, rmsd_flat, rmsf, rotation_mean, rotation_mean_with_options, superpose, superpose_with_options, symmetric, symmetric_with_options, torus_summary
from .ml import ArrowStream, AtomArrowTable, AtomTable, BondArrowTable, BondTable, ChainArrowTable, ChainTable, DLDataType, DLDevice, DLManagedTensor, DLTensor, Dataset, DatasetEntry, DatasetError, DatasetFilter, DatasetSplit, DatasetWarning, DlpackError, DlpackTensor, ExportCost, GraphError, LoadError, ManifestEntry, PdbioxExtension, ResidueArrowTable, ResidueTable, SplitOptions, SplitRatios, SplitStrategy, TableFileError, extension_name, graph, write_atom_ipc, write_atom_ipc_with_metadata, write_atom_parquet, write_atom_parquet_with_metadata
from .query import AltlocPolicy, AnalysisPolicy, AssemblyChoice, Evaluation, LogicalPlan, MissingPolicy, ModelChoice, Namespace, PhysicalQuery, Query, QueryBuilder, Selection, col
from .seq import *
from .compare import *
from .spatial import AtomsWithin, AutoBackendProfile, CellGridOptions, KdPeriodicOptions, NeighborList, NeighborListOptions, NeighborPair, NeighborPairs, NeighborSkinProfile, NeighborTable, PeriodicBox, PeriodicImage, SpatialError, SpatialOption, SpatialPlan, SpatialSearchOptions, atoms_within, atoms_within_with_options, nearest_neighbors, neighbor_pairs, neighbor_pairs_with_options
from .core.contract import Analysis, Assumption, AssumptionSource, Coverage, Diagnostic, ImpactEstimate, Provenance, Status
from .analysis import analyse_chain_interface, analyse_contacts, analyse_half_sphere_exposure, analyse_nucleic_torsions
from .surface import SurfaceWorkflowOptions, SurfaceWorkflowResult, analyse_surface_geometry
from .traj import CartesianFit, analyse_diffusion, analyse_pca, analyse_torsion_pca
from .validate import validate_bond_lengths, validate_cis_peptides, validate_clashes, validate_completeness, validate_planarity, validate_quality, validate_valence
from ._trajectory import DiffusionMap, DmsBond, DmsCell, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion, PcaResult, PeriodicAngle, Rotation3, SurfaceMesh, Trajectory, TrajectoryFormat, TrajectoryUnits, TrajectoryWriteOptions, read_dms, write_dms
from .modelcif import GlobalMetric, LocalMetric, MetricDefinition, ModelCategory, ModelCif, ModelDescription, ModelRow, PairwiseMetric, ProtocolStep, QualityMetrics, SoftwareGroup, Target, Template, lower_model_cif
from .pdb import MMTF_METADATA_EXTENSION, MmtfEntityMetadata, MmtfGroupMetadata, MmtfMetadata, MmtfOptionalField, PDB_HEADERS_EXTENSION, PdbHeaderRecord, PdbHeaders, PdbIdentifierNamespace, PdbOptions, PdbWriteOptions, mmtf_metadata, with_mmtf_metadata, write_mmtf_with_metadata

class Plan:
    def __init__(self, **operations: object) -> None: ...
    def execute(self, structure: Structure | None = None) -> dict[str, object]: ...
    def execute_report(self, structure: Structure | None = None) -> PlanResult: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> Plan: ...
    def __repr__(self) -> str: ...

@final
class PlanResult:
    results: dict[str, object]
    cached_index_count: int

@final
class Contacts:
    def __init__(self, *, left: str, right: str, cutoff: float = ..., backend: SpatialBackend | None = None, policy: AnalysisPolicy | None = ...) -> None: ...
    def execute(self, structure: Structure) -> Analysis: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> Contacts: ...
    def __repr__(self) -> str: ...

@final
class Rmsd:
    def __init__(self, *, mobile: NDArray[float32], reference: NDArray[float32]) -> None: ...
    def execute(self) -> float: ...
    def explain(self) -> dict[str, object]: ...
    def to_dict(self) -> dict[str, object]: ...
    @staticmethod
    def from_dict(config: dict[str, object]) -> Rmsd: ...
    def __repr__(self) -> str: ...
@final
class Limits:
    def __init__(self, decompressed_bytes: int, compression_ratio: int, rows_per_category: int, nesting_depth: int, dictionary_entries: int) -> None: ...
    @staticmethod
    def standard() -> Limits: ...
@final
class ReadScope:
    def __init__(self, only_first_model: bool, only_atomic_coords: bool, discard_hydrogens: bool) -> None: ...
    @staticmethod
    def all() -> ReadScope: ...
@final
class ReadOptions:
    def __init__(self, format: Format, mode: ParseMode, scope: ReadScope, missing_element_policy: MissingElementPolicy, ambiguous_residue_boundary_policy: AmbiguousResidueBoundaryPolicy, limits: Limits) -> None: ...
    @staticmethod
    def standard() -> ReadOptions: ...
@final
class PolicyOverrides:
    def __init__(self, *, assembly: str | None = ..., model: str | None = ..., altloc: str | None = ..., identifiers: str | None = ..., missing_atoms: str | None = ..., hydrogens: str | None = ..., atom_equivalence: str | None = ..., symmetry: str | None = ..., alignment: str | None = ..., precision: str | None = ..., periodic: str | None = ..., vdw_radii: str | None = ..., contact_def: str | None = ..., float_tolerance_relative: float | None = ..., float_tolerance_absolute: float | None = ...) -> None: ...
    assembly: str | None
    model: str | None
    altloc: str | None
    identifiers: str | None
    missing_atoms: str | None
    hydrogens: str | None
    atom_equivalence: str | None
    symmetry: str | None
    alignment: str | None
    precision: str | None
    periodic: str | None
    vdw_radii: str | None
    contact_def: str | None
    float_tolerance_relative: float | None
    float_tolerance_absolute: float | None
    def apply_to(self, policy: AnalysisPolicy | None = None) -> AnalysisPolicy: ...
@final
class OutputConfiguration:
    def __init__(self, format: str | None = ...) -> None: ...
    format: str | None
@final
class ChemistryConfiguration:
    def __init__(self, ccd_cache: PathLike[str] | None = ..., ccd_version: str | None = ...) -> None: ...
    ccd_cache: PathLike[str] | None
    ccd_version: str | None
@final
class ApplicationConfiguration:
    policy: PolicyOverrides
    output: OutputConfiguration
    chem: ChemistryConfiguration
@final
class ReadReport:
    structure: Structure
    findings: list[str]
@final
class CountDifference:
    left: int
    right: int
@final
class MetadataDifference:
    id: tuple[str | None, str | None] | None
    title: tuple[str | None, str | None] | None
    method: tuple[str | None, str | None] | None
    resolution: tuple[float | None, float | None] | None
    cell: tuple[tuple[list[float], list[float]] | None, tuple[list[float], list[float]] | None] | None
@final
class StructureDifference:
    metadata: MetadataDifference
    models: CountDifference | None
    entities: CountDifference | None
    chains: CountDifference | None
    residues: CountDifference | None
    atoms: CountDifference | None
    bonds: CountDifference | None
    changed_entities: int
    changed_chains: int
    changed_residues: int
    changed_atoms: int
    changed_bonds: int
    changed_positions: int
    maximum_displacement: float | None
    def is_empty(self) -> bool: ...
def structure_difference(left: Structure, right: Structure, coordinate_tolerance: float = ..., options: StructureDifferenceOptions | None = ...) -> StructureDifference: ...
@final
class BondOrder:
    Single: BondOrder
    Double: BondOrder
    Triple: BondOrder
    Quadruple: BondOrder
    Aromatic: BondOrder
    Polymeric: BondOrder
    Unknown: BondOrder
@final
class BondProvenance:
    File: BondProvenance
    ChemicalComponentDictionary: BondProvenance
    InferredDistance: BondProvenance
    User: BondProvenance
@final
class BondRecord:
    def __init__(self, atom_a: AtomIndex, atom_b: AtomIndex, order: BondOrder = ..., provenance: BondProvenance = ...) -> None: ...
    atom_a: int
    atom_b: int
    order: BondOrder
    provenance: BondProvenance
    atom_a_index: AtomIndex
    atom_b_index: AtomIndex
@final
class BondInference:
    def __init__(self, scale: float = 1.15, lower_bound: float = 0.4, *, exclude_across_chains: bool = False, respect_existing: bool = True, backend: SpatialBackend = SpatialBackend.Auto) -> None: ...
    @staticmethod
    def standard() -> BondInference: ...
    scale: float
    lower_bound: float
    exclude_across_chains: bool
    respect_existing: bool
    backend: SpatialBackend
@final
class BondInferenceReport:
    structure: Structure
    skipped_atoms: list[int]

@final
class ProteinAlphaTrace:
    chain: int
    positions: list[_Coordinate | None]

@final
class BackboneTorsionRecord:
    residue: int
    torsions: BackboneTorsions

@final
class SideChainTorsionRecord:
    residue: int
    atoms: list[str]
    torsions: list[float | None]

@final
class SideChainTorsionReport:
    records: list[SideChainTorsionRecord]
    findings: list[str]
    dictionary_version: str
@final
class Hedron:
    @staticmethod
    def from_points(i: _Coordinate, j: _Coordinate, k: _Coordinate) -> Hedron | None: ...
    first_length: float
    angle: float
    second_length: float
@final
class Dihedron:
    @staticmethod
    def from_points(i: _Coordinate, j: _Coordinate, k: _Coordinate, l: _Coordinate) -> Dihedron | None: ...
    first: Hedron
    angle: float
    length: float
    torsion: float
@final
class InternalAtom:
    atom: int
    references: tuple[int, int, int]
    coordinate: Dihedron
@final
class BatFrame:
    seed_positions: list[_Coordinate]
    coordinates: list[tuple[float, float, float]]
@final
class InternalCoordinates:
    seeds: list[tuple[int, _Coordinate]]
    atoms: list[InternalAtom]
    def measure_bat(self, positions: list[_Coordinate | None]) -> BatFrame: ...
    def rebuild_bat(self, frame: BatFrame) -> list[_Coordinate | None]: ...
    def rebuild(self) -> list[_Coordinate | None]: ...
@final
class CifValue:
    kind: str
    text: str | None
    integer: int | None
    number: float | None
    def as_identifier(self) -> str | None: ...
    def is_recorded(self) -> bool: ...
@final
class Column:
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> CifValue | None: ...
    values: list[CifValue]
@final
class Category:
    name: str
    row_count: int
    items: list[str]
    def column(self, item: str) -> Column | None: ...
    def value(self, item: str, row: int) -> CifValue | None: ...
    def text(self, item: str, row: int) -> str | None: ...
    def identifier(self, item: str, row: int) -> str | None: ...
@final
class DataBlock:
    name: str
    categories: list[Category]
    def category(self, name: str) -> Category | None: ...
@final
class Document:
    def __len__(self) -> int: ...
    is_empty: bool
    blocks: list[DataBlock]
    def first_block(self) -> DataBlock | None: ...
@final
class Model:
    @property
    def index(self) -> int: ...
    @property
    def number(self) -> int | None: ...
    @property
    def index_id(self) -> ModelIndex: ...
    @property
    def chains(self) -> Chains: ...
    def chain_at(self, position: int) -> Chain: ...
    def chain(self, label: str) -> Chain: ...
@final
class Models:
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Model: ...
@final
class Chain:
    @property
    def index(self) -> int: ...
    @property
    def index_id(self) -> ChainIndex: ...
    @property
    def label(self) -> str | None: ...
    @property
    def label_asym_id(self) -> SymbolId | None: ...
    @property
    def auth_label(self) -> str | None: ...
    @property
    def auth_asym_id(self) -> SymbolId | None: ...
    @property
    def residues(self) -> Residues: ...
    def residue_at(self, position: int) -> Residue: ...
    def residue(self, number: int) -> Residue | None: ...
    def entity(self) -> int | None: ...
    def entity_kind(self) -> EntityKind: ...
    def polymer_kind(self) -> PolymerKind: ...
    def observed_sequence(self) -> list[str]: ...
    def canonical_sequence(self) -> list[str]: ...
    def reference_sequences(self) -> list[ReferenceSequence]: ...
    def sequence_mapping(self) -> list[SequenceMapping]: ...
    def missing_residues(self) -> list[MissingResidue]: ...
@final
class Chains:
    def __len__(self) -> int: ...
    def __getitem__(self, key: int | str) -> Chain: ...
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...
@final
class Residue:
    @property
    def index(self) -> int: ...
    @property
    def index_id(self) -> ResidueIndex: ...
    @property
    def name(self) -> str | None: ...
    @property
    def label_comp_id(self) -> SymbolId | None: ...
    @property
    def auth_comp_id(self) -> SymbolId | None: ...
    @property
    def auth_name(self) -> str | None: ...
    @property
    def label_seq_id(self) -> int | None: ...
    @property
    def auth_seq_id(self) -> int | None: ...
    @property
    def ins_code(self) -> str | None: ...
    @property
    def is_het(self) -> bool: ...
    @property
    def atoms(self) -> ResidueAtoms: ...
    def atom_at(self, position: int) -> Atom: ...
    def atom(self, name: str) -> Atom: ...
@final
class Residues:
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Residue: ...
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...
@final
class ResidueAtoms:
    def __len__(self) -> int: ...
    def __getitem__(self, key: int | str) -> Atom: ...
@final
class Atom:
    @property
    def index(self) -> int: ...
    @property
    def index_id(self) -> AtomIndex: ...
    @property
    def name(self) -> str | None: ...
    @property
    def name_symbol(self) -> SymbolId | None: ...
    @property
    def auth_name(self) -> str | None: ...
    @property
    def auth_name_symbol(self) -> SymbolId | None: ...
    @property
    def altloc(self) -> str | None: ...
    def alt_label(self) -> str | None: ...
    @property
    def alt_id(self) -> AltId | None: ...
    @property
    def component(self) -> str | None: ...
    def component_name(self) -> str | None: ...
    @property
    def component_id(self) -> SymbolId | None: ...
    @property
    def element(self) -> str | None: ...
    @property
    def element_value(self) -> Element | None: ...
    @property
    def coord(self) -> _Coordinate | None: ...
    def position(self) -> _Coordinate | None: ...
    @property
    def b_factor(self) -> float | None: ...
    @property
    def occupancy(self) -> float | None: ...
    @property
    def formal_charge(self) -> int | None: ...
    @property
    def atom_site_id(self) -> int | None: ...
    @property
    def residue(self) -> Residue | None: ...

@final
class Atoms:
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Atom: ...
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...
@final
class Bonds:
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> BondRecord: ...
    @property
    def available(self) -> bool: ...
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...
@final
class CoordinateEdit:
    def __enter__(self) -> NDArray[float32]: ...
    def __exit__(self, exc_type, exc_value, traceback) -> bool: ...
