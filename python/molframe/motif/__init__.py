"""Declarative functional geometry: motifs, mappings and verdict profiles."""

from .._native import motif as _native

ConstraintResult = _native.ConstraintResult
MotifEvaluation = _native.MotifEvaluation
MotifReport = _native.MotifReport
evaluate = _native.evaluate
__all__ = ["ConstraintResult", "MotifEvaluation", "MotifReport", "evaluate"]
