from typing import final
from . import Format, ReadOptions, ReadReport, Structure, read_mmtf, read_pdb, read_pdbqt, read_pqr, write_mmtf, write_pdb, write_pdbqt, write_pqr
from . import fixed, hybrid36

@final
class PdbWriteOptions:
    def __init__(self, chain_map: dict[str, str] | None = ..., hybrid36: bool = ..., namespace: PdbIdentifierNamespace = ...) -> None: ...
@final
class PdbIdentifierNamespace:
    Label: PdbIdentifierNamespace
    Auth: PdbIdentifierNamespace
PdbOptions = PdbWriteOptions
PDB_HEADERS_EXTENSION: str
MMTF_METADATA_EXTENSION: str
@final
class MmtfOptionalField:
    BFactor: MmtfOptionalField
    Occupancy: MmtfOptionalField
    AtomId: MmtfOptionalField
    AltLoc: MmtfOptionalField
    InsCode: MmtfOptionalField
    SequenceIndex: MmtfOptionalField
    ChainName: MmtfOptionalField
    EntityList: MmtfOptionalField
@final
class MmtfGroupMetadata:
    def __init__(self, name: str, atom_names: list[str], elements: list[str] | None, single_letter_code: str, chem_comp_type: str) -> None: ...
    name: str
    atom_names: list[str]
    elements: list[str] | None
    single_letter_code: str
    chem_comp_type: str
@final
class MmtfEntityMetadata:
    def __init__(self, description: str, kind: str, sequence: str) -> None: ...
    description: str
    kind: str
    sequence: str
@final
class MmtfMetadata:
    def __init__(self, space_group: str | None = ..., groups: list[MmtfGroupMetadata] = ..., entities: list[MmtfEntityMetadata] = ..., optional_fields: list[MmtfOptionalField] = ...) -> None: ...
    space_group: str | None
    groups: list[MmtfGroupMetadata]
    entities: list[MmtfEntityMetadata]
    optional_fields: list[MmtfOptionalField]
@final
class PdbHeaderRecord:
    name: str
    line: str
@final
class PdbHeaders:
    def __len__(self) -> int: ...
    is_empty: bool
    records: list[PdbHeaderRecord]
    def named(self, name: str) -> list[PdbHeaderRecord]: ...
    classification: str | None
    deposition_date: str | None

def mmtf_metadata(structure: Structure) -> MmtfMetadata | None: ...
def with_mmtf_metadata(structure: Structure, metadata: MmtfMetadata) -> Structure: ...
def write_mmtf_with_metadata(structure: Structure, metadata: MmtfMetadata) -> bytes: ...
def read(data: bytes | str, options: ReadOptions | None = None) -> ReadReport: ...
def write(structure: Structure, options: PdbWriteOptions) -> str: ...
def write_selected(structure: Structure, options: PdbWriteOptions, selection: object | None = ...) -> str: ...
def pdb_write_selected(structure: Structure, options: PdbWriteOptions, selection: object | None = ...) -> str: ...

def field_text(line: str, start: int, stop: int) -> str: ...
def field_raw(line: str, start: int, stop: int) -> str: ...
def field_record(line: str) -> str: ...
def field_integer(line: str, start: int, stop: int) -> int | None: ...
def field_real(line: str, start: int, stop: int) -> float | None: ...
def hybrid36_decode(field: str, width: int) -> int | None: ...
def hybrid36_encode(value: int, width: int) -> str | None: ...
def hybrid36_needs_encoding(value: int, width: int) -> bool: ...
