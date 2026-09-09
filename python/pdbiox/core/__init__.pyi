from typing import Any
from .._io_types import CapacityError, TableError

from .contract import Analysis, Assumption, AssumptionSource, Coverage, Diagnostic, ImpactEstimate, Provenance, Status
from ..chem import Element
from .metadata import EntityKind, PolymerKind
from . import (
    annotation, bond, chunk, column, contract, coords, diagnostic, element,
    execution, index, io, limits, optional, parallel, provider, selection, span, structure, symbol, topology,
)
from ..query import AnalysisPolicy, Selection

from .._native import (
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION,
    COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION, HBOND_ACCEPTOR_ANNOTATION,
    HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION, PLDDT_ANNOTATION,
    POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION, STEREO_CONFIGURATION_ANNOTATION,
    Aabb, AltId, AmbiguousResidueBoundaryPolicy, AnnotationColumn, AtomAnnotation, AtomAnnotations, AtomChunk, AtomChunkStats,
    ChunkId, CountDifference, MetadataDifference, StructureDifference, structure_difference,
    AtomIndex, AtomRecord, BitVec, ChainRecord, ChainTable, ChunkBuilder, ColumnKind, CoordinateStore,
    ElementMask, EncodedColumn, EntityTable, ExtensionStore, Extremes,
    BondAdjacency, BondIndex, BondOrder, BondProvenance,
    BondRecord, BondTable, BondTableBuilder, ByteSpan, ChainIndex, Compression, CoordinateBlock,
    OutputOptions, DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
    Format, InputBuffer, InputKind, Limits, MissingElementPolicy, ParseMode, ReadOptions, ReadReport,
    ReadResult, ReadScope, Reader, Select, SelectAll,
    CoordinateGeneration, DictionaryFull, EntityIndex, InstanceId, Interner, ModelIndex, ModelTable,
    OptionalI32, OptionalSymbol, ParentMapping, Position, Presence, ResidueIndex, MissingResidue, ResidueRecord,
    ResidueTable, Structure, StructureData, StructureView, StructureEditor, CoordinateEditor,
    StructureBatch, StructureBatchReader, collect_structure, open_structure_batches,
    Class, Code, ContextItem, Diagnostics, Kind, Rendered, Severity, Strictness,
    SymbolId, TARGET_CHUNK_ATOMS, Topology, ValidityMask,
    bit_width, pack, write_output,
    unpack_one,
)

from .._core_primitives import *
from .._provider import *
from .._core_topology import *
from .parallel import BlockPlan, ReductionPolicy, DEFAULT_BLOCK_ITEMS
from .execution import (
    DEFAULT_MEMORY_BUDGET_BYTES, Backpressure,
    BatchDemand, CancellationToken, ExecutionContext, WindowedFile, MemoryBudget, MemoryLease,
    SpillArtifact, SpillFile, SpillReader,
    ScratchPolicy, TempStoragePolicy,
)

from .parallel import try_for_each_block_in
