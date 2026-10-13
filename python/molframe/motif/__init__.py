"""Declarative functional geometry: motifs, mappings and verdict profiles."""

from .._native import motif as _native

ConstraintResult = _native.ConstraintResult
MotifEvaluation = _native.MotifEvaluation
MotifReport = _native.MotifReport
MotifMatches = _native.MotifMatches
MatchIterator = _native.MatchIterator
evaluate = _native.evaluate
count = _native.count
matches = _native.matches
first = _native.first
iter_matches = _native.iter_matches
__all__ = [
    "ConstraintResult",
    "MatchIterator",
    "MotifEvaluation",
    "MotifMatches",
    "MotifReport",
    "count",
    "evaluate",
    "first",
    "iter_matches",
    "matches",
]
