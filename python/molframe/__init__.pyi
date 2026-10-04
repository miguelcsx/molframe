from collections.abc import Mapping, Sequence
from os import PathLike
from typing import ClassVar, Generic, Literal, Protocol, TypeVar

from . import (
    analysis as analysis,
    audit as audit,
    chemistry as chemistry,
    compare as compare,
    crystal as crystal,
    errors as errors,
    formats as formats,
    geometry as geometry,
    interop as interop,
    motif as motif,
    query as query,
    sel as sel,
    sequence as sequence,
    spatial as spatial,
    surface as surface,
    trajectory as trajectory,
    validation as validation,
)
from ._structure import (
    Annotation as Annotation,
    Annotations as Annotations,
    Atom as Atom,
    Atoms as Atoms,
    BondTable as BondTable,
    Chain as Chain,
    Chains as Chains,
    CoordinateEditor as CoordinateEditor,
    Entities as Entities,
    Entity as Entity,
    EntryMetadata as EntryMetadata,
    Model as Model,
    Models as Models,
    Query as Query,
    Residue as Residue,
    Residues as Residues,
    ResidueSelection as ResidueSelection,
    Selection as Selection,
    Structure as Structure,
    StructureEditor as StructureEditor,
)
from .errors import (
    Cancelled as Cancelled,
    ConsistencyError as ConsistencyError,
    ConversionError as ConversionError,
    Diagnostic as Diagnostic,
    GeometryError as GeometryError,
    IndeterminateError as IndeterminateError,
    InternalError as InternalError,
    MemoryBudgetError as MemoryBudgetError,
    MolframeError as MolframeError,
    MolframeIndexError as MolframeIndexError,
    MolframeIOError as MolframeIOError,
    MolframeKeyError as MolframeKeyError,
    MolframeTypeError as MolframeTypeError,
    MolframeValueError as MolframeValueError,
    MolframeWarning as MolframeWarning,
    ParseError as ParseError,
    PolicyError as PolicyError,
    QueryError as QueryError,
    QueryWarning as QueryWarning,
    ResourceError as ResourceError,
    SchemaError as SchemaError,
)

_T_co = TypeVar("_T_co", covariant=True)

class Float32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Float64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class UInt32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class UInt64Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class BoolArray(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Int32Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class UInt8Array(Protocol):
    @property
    def shape(self) -> tuple[int, ...]: ...

class Coverage:
    @property
    def intended(self) -> int: ...
    @property
    def used(self) -> int: ...
    @property
    def missing(self) -> int: ...
    @property
    def ambiguous(self) -> int: ...
    @property
    def fraction(self) -> float: ...

class Analysis(Generic[_T_co]):
    @property
    def value(self) -> _T_co:
        """The answer; raises ``IndeterminateError`` when the policy leaves none."""

    @property
    def is_determinate(self) -> bool:
        """True when there is an answer to read."""

    @property
    def indeterminacy(self) -> str | None:
        """Why there is no answer, or ``None`` when there is one."""

    @property
    def atom_origin(self) -> UInt32Array | None:
        """For each atom analysed, the input atom it is; ``None`` when they are the input's own.

        The atom indices in a result number the system the policy built (some atoms dropped,
        or copied for an assembly), and copies of an atom share its origin.
        """

    @property
    def policy_reads(self) -> list[str] | None:
        """The policy decisions the analysis applied, or ``None`` when it did not record them."""

    @property
    def status(self) -> Literal["complete", "partial", "ambiguous", "indeterminate"]: ...
    @property
    def coverage(self) -> Coverage: ...
    @property
    def warnings(self) -> list[str]: ...
    @property
    def assumptions(self) -> list[str]: ...
    @property
    def provenance(self) -> str: ...
    @property
    def profile(self) -> str | None: ...

class Table:
    def __len__(self) -> int: ...
    @property
    def names(self) -> list[str]: ...
    def __getitem__(self, name: str) -> Float32Array | Float64Array | UInt32Array | UInt64Array: ...
    def __contains__(self, name: str) -> bool: ...
    def to_dict(self) -> dict[str, Float32Array | Float64Array | UInt32Array | UInt64Array]: ...
    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...

class AnalysisPolicy:
    def __init__(
        self,
        *,
        identifiers: Literal["label", "auth", "explicit"] | None = None,
        altloc: str | None = None,
        missing_atoms: Literal["ignore", "report", "indeterminate", "fail"] | None = None,
        hydrogens: Literal["explicit_only", "include_inferred", "exclude"] | None = None,
        symmetry: Literal["none", "crystallographic", "biological_assembly"] | None = None,
        vdw_radii: Literal["bondi", "amber_united", "charmm", "alvarez"] | None = None,
        precision: Literal["f32", "f64"] | None = None,
        assembly: str | None = None,
        model: str | None = None,
        atom_equivalence: Literal["none", "ccd", "explicit"] | None = None,
        alignment: str | None = None,
        periodic: Literal["none", "pbc", "minimum_image"] | None = None,
        contact_def: str | None = None,
        float_tolerance: tuple[float, float] | None = None,
    ) -> None:
        """Name the decisions an analysis would otherwise make silently.

        ``assembly`` is ``"asymmetric_unit"``, ``"biological:<id>"`` or
        ``"crystal:<radius>"``; ``model`` is ``"first"``, ``"all"``, ``"ensemble"`` or
        ``"index:<n>"``; ``alignment`` is ``"none"``, ``"global"``, ``"local"`` or
        ``"explicit:<selection>"``; ``contact_def`` is ``"distance:<tolerance>"`` or
        ``"surface:<probe>"``. Unnamed decisions keep the default profile's choice.
        """
    @property
    def identifiers(self) -> str: ...
    @property
    def altloc(self) -> str: ...
    @property
    def missing_atoms(self) -> str: ...
    @property
    def hydrogens(self) -> str: ...
    @property
    def symmetry(self) -> str: ...
    @property
    def vdw_radii(self) -> str: ...
    @property
    def precision(self) -> str: ...
    @property
    def assembly(self) -> str: ...
    @property
    def model(self) -> str: ...
    @property
    def atom_equivalence(self) -> str: ...
    @property
    def alignment(self) -> str: ...
    @property
    def periodic(self) -> str: ...
    @property
    def contact_def(self) -> str: ...
    @property
    def float_tolerance(self) -> tuple[float, float]: ...
    @property
    def profile(self) -> str | None: ...
    @property
    def fingerprint(self) -> str: ...

class SecondarySource:
    Unassigned: ClassVar[SecondarySource]
    File: ClassVar[SecondarySource]
    Dssp: ClassVar[SecondarySource]
    CaOnly: ClassVar[SecondarySource]
    def __int__(self) -> int: ...

class SecondaryStructure:
    Unknown: ClassVar[SecondaryStructure]
    Coil: ClassVar[SecondaryStructure]
    AlphaHelix: ClassVar[SecondaryStructure]
    Strand: ClassVar[SecondaryStructure]
    Turn: ClassVar[SecondaryStructure]
    ThreeTenHelix: ClassVar[SecondaryStructure]
    PiHelix: ClassVar[SecondaryStructure]
    OtherHelix: ClassVar[SecondaryStructure]
    BetaBridge: ClassVar[SecondaryStructure]
    Bend: ClassVar[SecondaryStructure]
    PolyProline: ClassVar[SecondaryStructure]
    def __int__(self) -> int: ...
    def is_helix(self) -> bool: ...
    def is_strand(self) -> bool: ...
    def is_sheet_like(self) -> bool: ...

class ExecutionContext:
    """Limits and a cancellation switch for the operations run under them.

    Omitted limits keep the library's bounded defaults. A context describes
    limits and can govern many calls at once; ``cancel()`` stops whichever are
    running, and they raise :class:`Cancelled`. Ctrl-C reaches any operation that
    takes ``context=`` as cancellation and surfaces as ``KeyboardInterrupt``.
    """

    def __init__(
        self,
        *,
        workers: int | None = None,
        memory_budget: int | None = None,
        scratch_bytes: int = 0,
        temp_directory: str | PathLike[str] | None = None,
        temp_bytes: int = 0,
        image_limit: int | None = None,
    ) -> None: ...
    def cancel(self) -> None: ...
    @property
    def is_cancelled(self) -> bool: ...
    @property
    def workers(self) -> int | None: ...
    @property
    def image_limit(self) -> int | None:
        """The most candidate symmetry images a crystal search may examine, or ``None``."""

    @property
    def memory_budget(self) -> int | None: ...
    @property
    def scratch_bytes(self) -> int: ...
    @property
    def temp_directory(self) -> str | None: ...
    @property
    def temp_bytes(self) -> int: ...

class QueryAliases:
    def __init__(self) -> None: ...
    def define(self, name: str, query: Query) -> Query | None: ...
    def remove(self, name: str) -> Query | None: ...
    def get(self, name: str) -> Query | None: ...
    @property
    def names(self) -> list[str]: ...
    def resolve(self, query: Query) -> Query: ...
    def __len__(self) -> int: ...
    def __contains__(self, name: str) -> bool: ...

type _Format = Literal[
    "auto", "mmcif", "pdbml", "bcif", "mmtf", "pdb", "pqr", "pdbqt", "sdf", "mol2", "smallcif"
]

class ReadOptions:
    """Every decision a read makes, stated once and reusable across reads.

    Omitted arguments keep the library's choices. ``only_categories`` and
    ``skip_categories`` choose which optional mmCIF/BinaryCIF categories or PDB
    records are read (names the structure cannot be built without are always read);
    the ``max_*`` limits bound hostile or oversized input.
    """

    def __init__(
        self,
        *,
        format: _Format = "auto",
        mode: Literal["strict", "permissive", "recover"] = "permissive",
        first_model_only: bool = False,
        coordinates_only: bool = False,
        discard_hydrogens: bool = False,
        missing_element: Literal["preserve_unknown", "infer_from_atom_name"] = "preserve_unknown",
        ambiguous_residue_boundary: Literal["reject", "infer_from_file_order"] = "reject",
        only_categories: Sequence[str] | None = None,
        skip_categories: Sequence[str] | None = None,
        max_decompressed_bytes: int | None = None,
        max_compression_ratio: int | None = None,
        max_rows_per_category: int | None = None,
        max_nesting_depth: int | None = None,
        max_dictionary_entries: int | None = None,
    ) -> None: ...
    @property
    def format(self) -> str: ...
    @property
    def mode(self) -> str: ...
    @property
    def first_model_only(self) -> bool: ...
    @property
    def coordinates_only(self) -> bool: ...
    @property
    def discard_hydrogens(self) -> bool: ...
    @property
    def missing_element(self) -> str: ...
    @property
    def ambiguous_residue_boundary(self) -> str: ...
    @property
    def only_categories(self) -> list[str] | None: ...
    @property
    def skip_categories(self) -> list[str] | None: ...

class Reader:
    def __init__(self, data: bytes, *, name: str | None = ...) -> None: ...
    @property
    def byte_length(self) -> int: ...
    def read(
        self,
        *,
        format: _Format | None = None,
        options: ReadOptions | None = None,
        context: ExecutionContext | None = None,
    ) -> Structure: ...

class StructureBatch:
    """One bounded run of atom rows, in file order.

    ``occupancies`` and ``b_factors`` are ``NaN`` where the file records no value.
    """

    def __len__(self) -> int: ...
    @property
    def models(self) -> Int32Array: ...
    @property
    def sequences(self) -> Int32Array: ...
    @property
    def elements(self) -> UInt8Array: ...
    @property
    def positions(self) -> Float32Array: ...
    @property
    def occupancies(self) -> Float32Array: ...
    @property
    def b_factors(self) -> Float32Array: ...
    @property
    def atom_site_ids(self) -> UInt32Array: ...
    @property
    def heterogens(self) -> UInt8Array: ...
    @property
    def chains(self) -> list[str]: ...
    @property
    def components(self) -> list[str]: ...
    @property
    def atom_names(self) -> list[str]: ...
    @property
    def diagnostics(self) -> list[Diagnostic]: ...

class StructureBatches:
    def __iter__(self) -> StructureBatches: ...
    def __next__(self) -> StructureBatch: ...

class WorkflowNode: ...

class Workflow:
    def input(
        self,
        name: str,
        *,
        kind: Literal["coordinates", "float", "structure"],
        cost: Literal["borrow", "adopt", "decode", "copy", "materialize"] = ...,
    ) -> WorkflowNode: ...
    def centroid(self, coordinates: WorkflowNode) -> WorkflowNode: ...
    def rmsd(self, mobile: WorkflowNode, reference: WorkflowNode) -> WorkflowNode: ...
    def distance_matrix(self, coordinates: WorkflowNode) -> WorkflowNode: ...
    def atom_contacts(self, structure: WorkflowNode, cutoff: float) -> WorkflowNode: ...
    def output(self, name: str, node: WorkflowNode) -> None: ...
    def compile(self) -> CompiledWorkflow: ...

class CompiledWorkflow:
    def run(
        self,
        values: Mapping[str, object],
        *,
        copy: bool = ...,
        context: ExecutionContext | None = None,
    ) -> dict[str, object]: ...
    def explain(self) -> dict[str, object]: ...

__version__: str

def read(
    source: str | PathLike[str] | bytes,
    *,
    name: str | None = None,
    format: _Format | None = None,
    options: ReadOptions | None = None,
    context: ExecutionContext | None = None,
) -> Structure: ...
def read_with_diagnostics(
    source: str | PathLike[str] | bytes,
    *,
    name: str | None = None,
    format: _Format | None = None,
    options: ReadOptions | None = None,
    context: ExecutionContext | None = None,
) -> tuple[Structure, list[Diagnostic]]: ...
def open_structure_batches(
    path: str | PathLike[str],
    *,
    rows: int = 8192,
    bytes: int = 4_194_304,
    options: ReadOptions | None = None,
    context: ExecutionContext | None = None,
) -> StructureBatches: ...
