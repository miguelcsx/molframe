from typing import final

SEQUENCE_REFERENCES_EXTENSION: str

@final
class EntityKind:
    Polymer: EntityKind
    NonPolymer: EntityKind
    Water: EntityKind
    Branched: EntityKind
    Unknown: EntityKind
    def __repr__(self) -> str: ...

@final
class PolymerKind:
    @staticmethod
    def none() -> PolymerKind: ...
    Protein: PolymerKind
    Dna: PolymerKind
    Rna: PolymerKind
    NucleicHybrid: PolymerKind
    Saccharide: PolymerKind
    Other: PolymerKind
    def is_polymer(self) -> bool: ...
    def is_nucleic(self) -> bool: ...

@final
class EntryMetadata:
    id: str | None
    title: str | None
    method: str | None
    resolution: float | None

@final
class ReferenceSequence:
    id: str
    entity_id: str
    database_name: str | None
    database_code: str | None
    accession: str | None
    one_letter_code: str | None

@final
class ReferenceAlignment:
    id: str
    reference_id: str
    chain_ids: list[str]
    canonical: tuple[int, int]
    reference: tuple[int, int]

@final
class SequenceReferences:
    sequences: list[ReferenceSequence]
    alignments: list[ReferenceAlignment]

@final
class SequenceMapping:
    residue: int
    canonical_position: int | None
    reference_position: int | None
    reference_id: str | None
