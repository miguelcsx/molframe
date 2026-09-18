from typing import Any
from .._io_types import CapacityError, TableError

from .contract import Analysis, Assumption, AssumptionSource, Coverage, Diagnostic, ImpactEstimate, Provenance, Status
from ..chem import Element
from .metadata import EntityKind, EntryMetadata, PolymerKind, ReferenceAlignment, ReferenceSequence, SequenceMapping, SequenceReferences
from . import (
    annotation, bond, chunk, column, contract, coords, diagnostic, element,
    execution, index, io, limits, optional, parallel, provider, selection, span, structure, symbol, topology,
)
from ..query import AnalysisPolicy, Selection

from .._facade_models import (
    BondOrder, BondProvenance, BondRecord, CountDifference, Limits, MetadataDifference,
    ReadOptions, ReadScope, StructureDifference, structure_difference,
)

from .._facade_structure_io import (
    Structure, StructureBatch, StructureBatchReader, collect_structure,
    open_structure_batches,
)

from .._io_types import (
    AmbiguousResidueBoundaryPolicy, Format, MissingElementPolicy, ParseMode,
)

from .annotation import (
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION,
    COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION, HBOND_ACCEPTOR_ANNOTATION,
    HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION,
    PLDDT_ANNOTATION, POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
)

from .chunk import (
    TARGET_CHUNK_ATOMS,
)

from .._core_primitives import (
    Aabb, AnnotationColumn, AtomAnnotation, AtomAnnotations, AtomChunkStats, AtomSelection,
    BondAdjacency, BondTable, BondTableBuilder, ByteSpan, Class, Code,
    ColumnKind, Compression, ContextItem, CoordinateGeneration, DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES, Diagnostics,
    ElementMask, EncodedColumn, Extremes, InputBuffer, InputKind, Kind,
    OutputOptions, Position, ReadReport, ReadResult, Reader, Rendered,
    Select, SelectAll, Severity, Strictness, SymbolId,
)
from .._core_topology import (
    AltId, AtomChunk, AtomRecord, BitVec, ChainRecord, ChainTable,
    ChunkBuilder, CoordinateBlock, CoordinateEditor, CoordinateStore, DictionaryFull, EntityTable,
    ExtensionStore, Interner, MissingResidue, ModelTable, OptionalI32, OptionalSymbol,
    ParentMapping, Presence, ResidueRecord, ResidueTable, StructureData, StructureEditor,
    StructureView, Topology, ValidityMask, bit_width, pack, unpack_one,
)
from .._provider import (
    AtomEndpoint, BondChunk, BondChunkProvider, BondChunkRecord, ChunkDescriptor,
    ChunkId, ChunkLayout, DatasetCatalog, DatasetDescriptor, DatasetId, FrameChunk,
    FrameChunkProvider, LocalRow, LogicalRow, PayloadKind, PropertyChunk,
    PropertyChunkProvider, PropertyKind, PropertyValue, ProviderError, StructureChunk,
    StructureChunkProvider, TARGET_CHUNK_BONDS,
)
from .parallel import BlockPlan, ReductionPolicy, DEFAULT_BLOCK_ITEMS
from .execution import (
    DEFAULT_MEMORY_BUDGET_BYTES, Backpressure,
    BatchDemand, CancellationToken, ExecutionContext, WindowedFile, MemoryBudget, MemoryLease,
    SpillArtifact, SpillFile, SpillReader,
    ScratchPolicy, TempStoragePolicy,
)

from .parallel import try_for_each_block_in
