"""Type surface of the compiled extension module.

The extension registers its namespaces and classes; `molframe/__init__.py` and
each `molframe/<sub>/__init__.py` re-export them. Every name below *is* that
re-exported object, so the public stubs stay the single definition and this file
only points at them. Without it `molframe._native` resolves to an untyped
extension and strict mode reports `Unknown` for all ~200 re-exports.

Deliberately not a package: the runtime module is `molframe._native.abi3.so`, so
this is the sibling stub `molframe/_native.pyi`.
"""

from . import (
    Analysis as Analysis,
    AnalysisPolicy as AnalysisPolicy,
    Atom as Atom,
    Atoms as Atoms,
    BondTable as BondTable,
    Chain as Chain,
    Chains as Chains,
    CompiledWorkflow as CompiledWorkflow,
    Coverage as Coverage,
    ExecutionContext as ExecutionContext,
    Model as Model,
    Models as Models,
    Query as Query,
    QueryAliases as QueryAliases,
    Reader as Reader,
    ReadOptions as ReadOptions,
    Residue as Residue,
    Residues as Residues,
    ResidueSelection as ResidueSelection,
    SecondaryStructure as SecondaryStructure,
    Selection as Selection,
    Structure as Structure,
    StructureBatch as StructureBatch,
    StructureBatches as StructureBatches,
    StructureEditor as StructureEditor,
    Table as Table,
    Workflow as Workflow,
    WorkflowNode as WorkflowNode,
    __version__ as __version__,
    analysis as analysis,
    chemistry as chemistry,
    compare as compare,
    crystal as crystal,
    formats as formats,
    geometry as geometry,
    interop as interop,
    motif as motif,
    open_structure_batches as open_structure_batches,
    query as query,
    read as read,
    read_with_diagnostics as read_with_diagnostics,
    sel as sel,
    sequence as sequence,
    spatial as spatial,
    surface as surface,
    trajectory as trajectory,
    validation as validation,
)
