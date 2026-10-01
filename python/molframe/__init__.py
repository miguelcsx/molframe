
"""Curated Python contract for MolFrame."""

from . import _native
from ._native import (
    Analysis,
    AnalysisPolicy,
    Coverage,
    Atom,
    Atoms,
    Chain,
    Chains,
    CompiledWorkflow,
    Model,
    Models,
    Query,
    QueryAliases,
    QueryError,
    QueryWarning,
    Reader,
    Residue,
    ResidueSelection,
    Residues,
    Selection,
    Structure,
    StructureEditor,
    Table,
    Workflow,
    WorkflowNode,
    __version__,
    read,
)

geometry = _native.geometry
analysis = _native.analysis
trajectory = _native.trajectory
sequence = _native.sequence
crystal = _native.crystal
validation = _native.validation
motif = _native.motif
chemistry = _native.chemistry
compare = _native.compare
query = _native.query
sel = _native.sel
spatial = _native.spatial
surface = _native.surface
formats = _native.formats

__all__ = [
    "__version__",
    "Analysis",
    "AnalysisPolicy",
    "Coverage",
    "CompiledWorkflow",
    "Atom",
    "Atoms",
    "Chain",
    "Chains",
    "Model",
    "Models",
    "Query",
    "QueryAliases",
    "QueryError",
    "QueryWarning",
    "Reader",
    "Residue",
    "ResidueSelection",
    "Residues",
    "Selection",
    "Structure",
    "Table",
    "StructureEditor",
    "Workflow",
    "WorkflowNode",
    "analysis",
    "chemistry",
    "compare",
    "crystal",
    "formats",
    "geometry",
    "motif",
    "query",
    "read",
    "sel",
    "sequence",
    "spatial",
    "surface",
    "trajectory",
    "validation",
]
