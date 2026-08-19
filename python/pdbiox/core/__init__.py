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
    CapacityError, ChunkId, CountDifference, DictionaryVersion, DifferenceError, Fingerprint, MetadataDifference, ParameterValue, ProfileId, SourceRef,
    AtomIndex, AtomRecord, BitVec, ChainRecord, ChainTable, ChunkBuilder, ColumnKind, CoordinateStore,
    ElementMask, EntityTable, ExtensionStore, Extremes,
    BondAdjacency, BondIndex, BondOrder, BondProvenance,
    BondRecord, BondTable, BondTableBuilder, Compression, EncodedColumn, Format, InputBuffer, InputKind, Limits,
    MissingElementPolicy, ParseMode, ReadOptions, ReadReport, ReadResult, ReadScope, Reader, Select, SelectAll,
    ByteSpan, ChainIndex, CoordinateBlock, StructureDifference, StructureDifferenceOptions, TableError, ValueDifference, structure_difference,
    CoordinateGeneration, DictionaryFull, EntityIndex, InstanceId, Interner, ModelIndex, ModelTable,
    OptionalI32, OptionalSymbol, ParentMapping, Position, Presence, ResidueIndex, MissingResidue, Structure,
    SymbolId, TARGET_CHUNK_ATOMS, Topology, ResidueRecord, ResidueTable, StructureData,
    StructureView, StructureEditor, CoordinateEditor, ValidityMask, AtomSelection, Selection,
    Class, Code, ContextItem, Diagnostics, Kind, Rendered, Severity, Strictness,
    bit_width, pack, unpack_one, write_output,
)

__all__ = [
    "Analysis", "AnalysisPolicy", "Assumption", "AssumptionSource", "Coverage", "Diagnostic",
    "ImpactEstimate", "Provenance", "Status", "AlgorithmId", "AnalysisParameters",
    "CapacityError", "ChunkId", "CountDifference", "DictionaryVersion", "DifferenceError", "Fingerprint", "MetadataDifference", "ParameterValue", "ProfileId", "SourceRef",
    "AROMATIC_ATOM_ANNOTATION", "ATOM_RADIUS_ANNOTATION", "AUTODOCK_TYPE_ANNOTATION",
    "COMPONENT_KIND_ANNOTATION", "FORMAL_CHARGE_ANNOTATION", "HBOND_ACCEPTOR_ANNOTATION",
    "HBOND_DONOR_ANNOTATION", "PAE_ANNOTATION", "PARTIAL_CHARGE_ANNOTATION", "PLDDT_ANNOTATION",
    "POLYMER_ATOM_ROLE_ANNOTATION", "SEGMENT_ID_ANNOTATION", "STEREO_CONFIGURATION_ANNOTATION",
    "Aabb", "AltId", "AmbiguousResidueBoundaryPolicy", "AnnotationColumn", "AtomAnnotation", "AtomAnnotations", "AtomChunk",
    "Class", "Code", "ContextItem", "Diagnostics", "Kind", "Rendered", "Severity", "Strictness",
    "AtomChunkStats",
    "AtomIndex", "AtomRecord", "BitVec", "ChainRecord", "ChainTable", "ChunkBuilder",
    "ColumnKind", "EncodedColumn",
    "ElementMask", "EntityTable", "ExtensionStore", "Extremes",
    "BondAdjacency", "BondIndex", "BondOrder",
    "BondProvenance", "BondRecord", "BondTable", "BondTableBuilder", "ByteSpan", "ChainIndex",
    "Compression", "Format", "InputBuffer", "InputKind", "Limits", "MissingElementPolicy", "ParseMode",
    "ReadOptions", "ReadReport", "ReadResult", "ReadScope", "Reader",
    "StructureDifference", "StructureDifferenceOptions", "TableError", "ValueDifference", "structure_difference",
    "Select", "SelectAll",
    "CoordinateBlock", "CoordinateGeneration", "DictionaryFull", "EntityIndex", "InstanceId",
    "Interner", "ModelIndex", "ModelTable", "OptionalI32", "OptionalSymbol", "Position",
    "Presence", "ParentMapping", "ResidueIndex", "MissingResidue", "ResidueRecord", "ResidueTable", "Structure",
    "SymbolId", "TARGET_CHUNK_ATOMS", "Topology", "StructureData", "StructureView",
    "StructureEditor", "CoordinateEditor", "CoordinateStore", "AtomSelection", "Selection",
    "Element", "EntityKind", "EntryMetadata", "PolymerKind", "ReferenceAlignment", "ReferenceSequence", "SequenceMapping", "SequenceReferences",
    "ValidityMask", "bit_width", "pack", "unpack_one", "write_output",
    "annotation", "bond", "chunk", "column", "contract", "coords",
    "diagnostic", "element", "index", "io", "limits", "optional",
    "selection", "span", "structure", "symbol", "topology",
]

from . import (
    annotation, bond, chunk, column, contract, coords, diagnostic, element,
    index, io, limits, optional, selection, span, structure, symbol, topology,
)
