"""Core facade handles and immutable structure views."""

from .contract import Analysis, Assumption, AssumptionSource, Coverage, Diagnostic, ImpactEstimate, Provenance, Status
from ..chem import Element
from .metadata import EntityKind, EntryMetadata, PolymerKind, ReferenceAlignment, ReferenceSequence, SequenceMapping, SequenceReferences
from ..query import AnalysisPolicy
from .._native import (
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION,
    COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION, HBOND_ACCEPTOR_ANNOTATION,
    HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION, PLDDT_ANNOTATION,
    POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION, STEREO_CONFIGURATION_ANNOTATION,
    Aabb, AltId, AlgorithmId, AnalysisParameters, AmbiguousResidueBoundaryPolicy, AnnotationColumn, AtomAnnotation, AtomAnnotations, AtomChunk, AtomChunkStats,
    CapacityError, CountDifference, DictionaryVersion, DifferenceError, Fingerprint, MetadataDifference, ParameterValue, ProfileId, SourceRef,
    AtomIndex, AtomRecord, BitVec, ChainRecord, ChainTable, ChunkBuilder, ColumnKind, CoordinateStore,
    ElementMask, EntityTable, ExtensionStore, Extremes,
    BondAdjacency, BondIndex, BondOrder, BondProvenance,
    BondRecord, BondTable, BondTableBuilder, Compression, EncodedColumn, Format, InputBuffer, InputKind, Limits,
    OutputOptions, DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
    MemoryBudget, CancellationToken, ScratchPolicy, TempStoragePolicy,
    BlockPlan, ReductionPolicy, DEFAULT_BLOCK_ITEMS, try_for_each_block_in,
    ExecutionContext, WindowedFile, MemoryLease, SpillFile, SpillArtifact, SpillReader, BatchDemand, Backpressure,
    DatasetId, ChunkId, LogicalRow, LocalRow, PayloadKind, ChunkLayout,
    ChunkDescriptor, DatasetDescriptor, DatasetCatalog, PropertyKind, PropertyValue,
    AtomEndpoint, BondChunk, BondChunkRecord, StructureChunk, PropertyChunk, FrameChunk,
    BondChunkProvider, StructureChunkProvider,
    PropertyChunkProvider, FrameChunkProvider, ProviderError,
    DEFAULT_MEMORY_BUDGET_BYTES,
    MissingElementPolicy, ParseMode, ReadOptions, ReadReport, ReadResult, ReadScope, Reader, Select, SelectAll,
    ByteSpan, ChainIndex, CoordinateBlock, StructureDifference, StructureDifferenceOptions, TableError, ValueDifference, structure_difference,
    CoordinateGeneration, DictionaryFull, EntityIndex, InstanceId, Interner, ModelIndex, ModelTable,
    OptionalI32, OptionalSymbol, ParentMapping, Position, Presence, ResidueIndex, MissingResidue, Structure,
    SymbolId, TARGET_CHUNK_ATOMS, TARGET_CHUNK_BONDS, Topology, ResidueRecord, ResidueTable, StructureData,
    StructureView, StructureEditor, CoordinateEditor, ValidityMask, AtomSelection, Selection,
    StructureBatch, StructureBatchReader, collect_structure, open_structure_batches,
    Class, Code, ContextItem, Diagnostics, Kind, Rendered, Severity, Strictness,
    bit_width, pack, unpack_one, write_output,
)

from . import (
    annotation, bond, chunk, column, contract, coords, diagnostic, element,
    execution, index, io, limits, optional, parallel, provider, selection, span, structure, symbol, topology,
)

# The public surface is what this module binds, so it is read off the module
# rather than written out again: every name imported above is public, and a name
# that must stay private takes a leading underscore. The 42-line list this
# replaces said the same thing a second time, and nothing compared the two.
__all__ = sorted(name for name in dir() if not name.startswith("_"))
