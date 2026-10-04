"""Structure comparison operations."""

from .._native import compare as _native

CadScore = _native.CadScore
CeAlignment = _native.CeAlignment
ChainAssignment = _native.ChainAssignment
ContactSimilarity = _native.ContactSimilarity
DockQ = _native.DockQ
MappedComparison = _native.MappedComparison
assign_chains = _native.assign_chains
cad_score = _native.cad_score
ce_align = _native.ce_align
contact_areas = _native.contact_areas
contact_similarity = _native.contact_similarity
dockq = _native.dockq
gdt_ha = _native.gdt_ha
gdt_ts = _native.gdt_ts
lddt = _native.lddt
map_sequence_to_structure = _native.map_sequence_to_structure
mapped_dockq = _native.mapped_dockq
mapped_qs_score = _native.mapped_qs_score
qs_score = _native.qs_score
tm_score = _native.tm_score
weighted_rmsd = _native.weighted_rmsd

__all__ = [
    "CadScore",
    "CeAlignment",
    "ChainAssignment",
    "ContactSimilarity",
    "DockQ",
    "MappedComparison",
    "assign_chains",
    "cad_score",
    "ce_align",
    "contact_areas",
    "contact_similarity",
    "dockq",
    "gdt_ha",
    "gdt_ts",
    "lddt",
    "map_sequence_to_structure",
    "mapped_dockq",
    "mapped_qs_score",
    "qs_score",
    "tm_score",
    "weighted_rmsd",
]
