"""The import header the three facade stubs share."""

from os import PathLike

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

from .xtal import DEFAULT_CRYSTAL_IMAGE_LIMIT, DEFAULT_INSTANCE_LIMIT, SpaceGroup, SymmetryOperation, UnitCell, space_group_by_hall, space_group_by_number, space_group_by_symbol, space_group_setting, space_group_settings

from .xtal import ASSEMBLIES_EXTENSION, AffineTransform, AssemblyDef, AssemblyNeighbor, AssemblySet, AssemblyView, AtomInstance, CellTransform, ChainInstance, CrystalNeighbor, Generator, INSTANCE_ID_ANNOTATION, NCS_EXTENSION, NcsCode, NcsOperator, NcsSet, NcsView, OperExpression, Operator, Rational, SpaceGroupSetting, SymmetrySet, SYMMETRY_EXTENSION, collect_crystal_neighbors, lower_assemblies, lower_ncs, lower_symmetry

from .ml.graph import EdgeDirection, EdgeFeature, EdgeKind, Graph, GraphOptions, MissingFeaturePolicy, NodeFeature, NodeLevel, SpatialBackend

from .seq import *

from .compare import *

from .spatial import AtomsWithin, AutoBackendProfile, CellGridOptions, KdPeriodicOptions, NeighborList, NeighborListOptions, NeighborPair, NeighborPairs, NeighborSkinProfile, NeighborTable, PeriodicBox, PeriodicImage, SpatialError, SpatialOption, SpatialPlan, SpatialSearchOptions, atoms_within, atoms_within_with_options, nearest_neighbors, neighbor_pairs, neighbor_pairs_with_options

from .core.contract import Analysis, Assumption, AssumptionSource, Coverage, Diagnostic, ImpactEstimate, Provenance, Status

from .analysis import analyse_chain_interface, analyse_contacts, analyse_half_sphere_exposure, analyse_nucleic_torsions

from .surface import SurfaceWorkflowOptions, SurfaceWorkflowResult, analyse_surface_geometry

from .traj import CartesianFit, analyse_diffusion, analyse_pca, analyse_torsion_pca

from .validate import validate_bond_lengths, validate_cis_peptides, validate_clashes, validate_completeness, validate_planarity, validate_quality, validate_valence

from ._trajectory import DiffusionMap, DmsBond, DmsCell, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion, PcaResult, PeriodicAngle, Rotation3, SurfaceMesh, Trajectory, TrajectoryFormat, TrajectoryUnits, TrajectoryWriteOptions, write_dms

from .modelcif import GlobalMetric, LocalMetric, MetricDefinition, ModelCategory, ModelCif, ModelDescription, ModelRow, PairwiseMetric, ProtocolStep, QualityMetrics, SoftwareGroup, Target, Template, lower_model_cif

from .pdb import MMTF_METADATA_EXTENSION, MmtfEntityMetadata, MmtfGroupMetadata, MmtfMetadata, MmtfOptionalField, PDB_HEADERS_EXTENSION, PdbHeaderRecord, PdbHeaders, PdbIdentifierNamespace, PdbOptions, PdbWriteOptions, mmtf_metadata, with_mmtf_metadata, write_mmtf_with_metadata

