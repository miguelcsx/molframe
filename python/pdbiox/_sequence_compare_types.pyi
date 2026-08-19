"""Sequence, alignment, comparison, and mapping value types."""

from typing import final, Protocol
from numpy import float32, float64, uint32
from numpy.typing import NDArray
from . import Structure
from .core.contract import Analysis
from .chem import ComponentDictionary
from .geom import Rigid
from .query import AnalysisPolicy, Namespace

class Alphabet(Protocol):
    def symbols(self) -> bytes: ...
    def encode_symbol(self, symbol: int) -> int | None: ...
    def decode_symbol(self, code: int) -> int | None: ...
    def encode(self, symbols: bytes) -> bytes | None: ...
    def decode(self, codes: bytes) -> bytes | None: ...
    def len(self) -> int: ...
    def is_empty(self) -> bool: ...

Code = int

class AlphabetError(ValueError): ...
class SequenceError(ValueError): ...
class A2mError(ValueError): ...
class AlignError(ValueError): ...
class FastqError(ValueError): ...
class MatrixError(ValueError): ...
class MsaError(ValueError): ...
class RegionAlignError(ValueError): ...
class RegionError(ValueError): ...
class SequenceFormatError(ValueError): ...

class ProteinAlphabet:
    def __init__(self) -> None: ...
    def symbols(self) -> bytes: ...
    def encode_symbol(self, symbol: int) -> int | None: ...
    def decode_symbol(self, code: int) -> int | None: ...
    def encode(self, symbols: bytes) -> bytes | None: ...
    def decode(self, codes: bytes) -> bytes | None: ...
    def len(self) -> int: ...
    def is_empty(self) -> bool: ...

@final
class DnaAlphabet(ProteinAlphabet): ...
@final
class RnaAlphabet(ProteinAlphabet): ...
@final
class NucleotideAlphabet(ProteinAlphabet): ...

@final
class CustomAlphabet:
    def __init__(self, symbols: bytes) -> None: ...
    def symbols(self) -> bytes: ...
    def encode_symbol(self, symbol: int) -> int | None: ...
    def decode_symbol(self, code: int) -> int | None: ...
    def encode(self, symbols: bytes) -> bytes | None: ...
    def decode(self, codes: bytes) -> bytes | None: ...
    def len(self) -> int: ...
    def is_empty(self) -> bool: ...

@final
class Sequence:
    def __init__(self, alphabet: Alphabet, symbols: bytes) -> None: ...
    @staticmethod
    def from_codes(alphabet: Alphabet, codes: bytes) -> Sequence: ...
    alphabet: Alphabet
    codes: bytes
    symbols: bytes
    def len(self) -> int: ...
    def is_empty(self) -> bool: ...
    def __repr__(self) -> str: ...

PROTEIN: ProteinAlphabet
DNA: DnaAlphabet
RNA: RnaAlphabet
NUCLEOTIDE: NucleotideAlphabet

class KmerTableError(ValueError): ...
class SeedPatternError(ValueError): ...
class SimilarKmerError(ValueError): ...

@final
class KmerStorage:
    @staticmethod
    def exact() -> KmerStorage: ...
    @staticmethod
    def bucketed(buckets: int) -> KmerStorage: ...
    buckets: int | None
    is_exact: bool

@final
class SeedPattern:
    def __init__(self, mask: list[bool]) -> None: ...
    span: int
    weight: int

@final
class KmerTableOptions:
    def __init__(self, k: int, *, pattern: SeedPattern | None = ..., storage: KmerStorage | None = ...) -> None: ...
    k: int
    pattern: SeedPattern | None
    storage: KmerStorage

@final
class KmerHit:
    sequence: int
    position: int

@final
class KmerTable:
    @staticmethod
    def build(sequences: list[bytes], options: KmerTableOptions) -> KmerTable: ...
    def query(self, word: bytes) -> list[KmerHit]: ...
    word_length: int

@final
class Scoring:
    def __init__(self, match_score: int = 1, mismatch_score: int = -1, gap_open: int = -2, gap_extend: int = -1) -> None: ...
    @staticmethod
    def simple() -> Scoring: ...
    match_score: int
    mismatch_score: int
    gap_open: int
    gap_extend: int
    def substitution(self, left: int, right: int) -> int: ...

@final
class Alignment:
    @property
    def score(self) -> int: ...
    @property
    def columns(self) -> list[Column]: ...

@final
class AlignmentColumn:
    left: int | None
    right: int | None

Column = AlignmentColumn

@final
class MsaOptions:
    def __init__(self, scoring: Scoring, *, refinement_passes: int = 0) -> None: ...
    @staticmethod
    def progressive(scoring: Scoring) -> MsaOptions: ...
    def with_refinement_passes(self, refinement_passes: int) -> MsaOptions: ...

@final
class EmptyLddtPolicy:
    Perfect: EmptyLddtPolicy
    Error: EmptyLddtPolicy

@final
class LddtOptions:
    def __init__(self, inclusion_radius: float, minimum_reference_distance: float, tolerances: list[float], empty_policy: EmptyLddtPolicy) -> None: ...
    @staticmethod
    def standard(inclusion_radius: float) -> LddtOptions: ...

@final
class CeSignificanceProfile:
    OriginalWindowEight: CeSignificanceProfile

@final
class CeOptions:
    def __init__(self, window_size: int, max_gap: int, max_paths: int, fragment_similarity_threshold: float, path_similarity_threshold: float, significance: CeSignificanceProfile | None) -> None: ...
    @staticmethod
    def original() -> CeOptions: ...

@final
class CeAlignment:
    reference_indices: list[int]
    mobile_indices: list[int]
    fragment_count: int
    similarity: float
    z_score: float | None
    rmsd: float

@final
class EmptyQsPolicy:
    Perfect: EmptyQsPolicy
    Error: EmptyQsPolicy

@final
class QsOptions:
    def __init__(self, contact_distance: float, empty_policy: EmptyQsPolicy) -> None: ...
    @staticmethod
    def standard(contact_distance: float) -> QsOptions: ...

@final
class DockQOptions:
    def __init__(self, contact_distance: float, ligand_scale: float, interface_scale: float) -> None: ...

@final
class DockQ:
    fnat: float
    ligand_rmsd: float
    interface_rmsd: float
    score: float

@final
class ContactArea:
    def __init__(self, first: int, second: int, area: float) -> None: ...
    first: int
    second: int
    area: float

@final
class CadContact:
    first: int
    second: int
    reference_area: float
    model_area: float
    lost_area: float

@final
class LocalCad:
    residue: int
    reference_area: float
    lost_area: float
    score: float

@final
class CadScore:
    score: float
    reference_area: float
    lost_area: float
    contacts: list[CadContact]
    local: list[LocalCad]

@final
class ContactSimilarity:
    shared: int
    union: int
    jaccard: float

@final
class EquivalentAtomMapping:
    reference_to_model: list[int]

@final
class LigandRmsd:
    rmsd: float
    mapping: EquivalentAtomMapping

@final
class ChainSequence:
    label: str
    sequence: bytes

@final
class ChainMapping:
    reference: str
    target: str
    identity: float

@final
class ChainAlternative:
    reference: str
    target: str
    identity: float

@final
class ChainAssignment:
    primary: list[ChainMapping]
    alternatives: list[ChainAlternative]

@final
class ResidueMatch:
    query_position: int
    residue: int

@final
class FastaRecord:
    def __init__(self, id: str, description: str, sequence: bytes) -> None: ...
    id: str
    description: str
    sequence: bytes

@final
class FastqRecord:
    def __init__(self, id: str, description: str, sequence: bytes, quality: bytes) -> None: ...
    id: str
    description: str
    sequence: bytes
    quality: bytes

@final
class SequenceFormat:
    Fasta: SequenceFormat
    Fastq: SequenceFormat
    Stockholm: SequenceFormat
    Clustal: SequenceFormat
    Phylip: SequenceFormat
    A2m: SequenceFormat
    A3m: SequenceFormat
    Newick: SequenceFormat

@final
class SequenceDocumentKind:
    Records: SequenceDocumentKind
    Fastq: SequenceDocumentKind
    Tree: SequenceDocumentKind

@final
class SequenceDocument:
    @staticmethod
    def records(records: list[FastaRecord]) -> SequenceDocument: ...
    @staticmethod
    def fastq(records: list[FastqRecord]) -> SequenceDocument: ...
    @staticmethod
    def tree(tree: Tree) -> SequenceDocument: ...
    @property
    def kind(self) -> SequenceDocumentKind: ...
    def record_values(self) -> list[FastaRecord]: ...
    def fastq_values(self) -> list[FastqRecord]: ...
    def tree_value(self) -> Tree: ...

@final
class MatrixProfile:
    @staticmethod
    def blosum(level: int) -> MatrixProfile: ...
    @staticmethod
    def pam(distance: int) -> MatrixProfile: ...
    @staticmethod
    def identity() -> MatrixProfile: ...
    @staticmethod
    def nuc44() -> MatrixProfile: ...

@final
class SimilarKmerOptions:
    def __init__(self, minimum_score: int, max_results: int, max_candidates: int) -> None: ...

@final
class SimilarKmer:
    word: bytes
    score: int

@final
class AlignmentMode:
    Global: AlignmentMode
    Local: AlignmentMode
    SemiGlobal: AlignmentMode

@final
class PointMatch:
    def __init__(self, reference: int, model: int) -> None: ...
    reference: int
    model: int

@final
class PointMapping:
    def __init__(self, matches: list[PointMatch], reference_len: int, model_len: int) -> None: ...
    matches: list[PointMatch]
    def __len__(self) -> int: ...

@final
class ComparisonAlignment:
    @staticmethod
    def not_required() -> ComparisonAlignment: ...
    @staticmethod
    def rigid(transform: Rigid) -> ComparisonAlignment: ...
    is_required: bool
    transform: Rigid | None

@final
class DistanceMeasurement:
    def __init__(self, distances: list[float], rmsd: float) -> None: ...
    distances: list[float]
    rmsd: float

@final
class ComparisonVerdict:
    def __init__(self, maximum_rmsd: float, passed: bool) -> None: ...
    maximum_rmsd: float
    passed: bool

@final
class InterfaceRmsd:
    def __init__(self, value: float) -> None: ...
    value: float

@final
class PocketRmsd:
    def __init__(self, value: float) -> None: ...
    value: float

@final
class RegionOptions:
    def __init__(self, left: tuple[int, int], right: tuple[int, int], mode: AlignmentMode, *, band: int | None = None) -> None: ...

@final
class SubstitutionMatrix:
    def score(self, left: int, right: int) -> int: ...
    def get(self, left: int, right: int) -> int: ...
    def identity(self) -> MatrixIdentity: ...
    name: str
    version: str
    source: str
    retrieved: str

@final
class MatrixIdentity:
    name: str
    version: str
    source: str
    retrieved: str
