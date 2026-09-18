"""GROMACS TPR topology format."""

from .._trajectory import TprAtom, TprBond, TprError, TprHeader, TprResidue, TprTopology, parse_tpr

__all__ = ["TprAtom", "TprBond", "TprError", "TprHeader", "TprResidue", "TprTopology", "parse_tpr"]
