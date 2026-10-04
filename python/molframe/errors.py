"""The exceptions MolFrame raises, and the findings they carry.

Every failure is a :class:`MolframeError`. The subclass says what kind of
failure it was, and where Python has an obvious built-in for it the exception
is also an instance of that built-in, so ``except ValueError`` keeps working
for a bad value and ``except MemoryError`` for an exhausted budget.

Each exception exposes the structured fields of the diagnostic behind it:
``code`` (such as ``"MOLFRAME-E5101"``), ``message``, ``remedy`` (what to do
about it), ``span`` (a byte range in the input, when there is one) and
``findings`` (every diagnostic of a failed read or query, in order).
Non-fatal diagnostics are :class:`MolframeWarning` instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence


@dataclass(frozen=True, slots=True)
class Diagnostic:
    """One finding about data: a registered code, what happened and what to do."""

    code: str
    message: str
    remedy: str | None = None
    span: tuple[int, int] | None = None
    severity: str | None = None
    context: dict[str, str] = field(default_factory=dict[str, str])


class MolframeError(Exception):
    """Base class of every error MolFrame raises."""

    code: str | None
    message: str
    remedy: str | None
    span: tuple[int, int] | None
    findings: tuple[Diagnostic, ...]

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        remedy: str | None = None,
        span: tuple[int, int] | None = None,
        findings: Sequence[Diagnostic] = (),
    ) -> None:
        """Record the message and the structured fields of the diagnostic."""
        super().__init__(message)
        self.message = message
        self.code = code
        self.remedy = remedy
        self.span = span
        self.findings = tuple(findings)

    def __str__(self) -> str:
        """Return the message and its code, without the quoting a bare ``KeyError`` adds."""
        if self.code is None:
            return self.message
        detail = ""
        if self.findings and self.findings[0].context:
            detail = " (" + ", ".join(f"{k}={v}" for k, v in self.findings[0].context.items()) + ")"
        return f"{self.message}{detail} [{self.code}]"


class MolframeValueError(MolframeError, ValueError):
    """An argument or input value is not acceptable."""


class MolframeKeyError(MolframeError, KeyError):
    """A name or label that was asked for does not exist."""


class MolframeIndexError(MolframeError, IndexError):
    """A position that was asked for is out of range."""


class MolframeTypeError(MolframeError, TypeError):
    """An argument has the wrong type."""


class ParseError(MolframeValueError):
    """The text or bytes of a file are malformed (diagnostic class 1xxx)."""


class SchemaError(MolframeValueError):
    """A file does not conform to its schema or dictionary (class 2xxx)."""


class ConsistencyError(MolframeValueError):
    """A structure is internally inconsistent (class 3xxx)."""


class ConversionError(MolframeValueError):
    """Data cannot be converted or selected as asked (class 4xxx)."""


class QueryError(ConversionError):
    """A selection query that cannot be compiled or evaluated."""


class GeometryError(MolframeValueError):
    """Geometry or numerics refused the input (class 5xxx)."""


class PolicyError(MolframeValueError):
    """The analysis policy or request is contradictory or unsupported (class 6xxx)."""


class IndeterminateError(MolframeError):
    """An analysis has no defensible answer under its policy, and its value was asked for.

    An indeterminate :class:`Analysis` holds a reason instead of a value, so the way
    to a number is ``analysis.value`` and it raises this. Check ``analysis.is_determinate``
    first, or read ``analysis.indeterminacy`` for why.
    """


class ResourceError(MolframeError, RuntimeError):
    """A resource limit, input or output failure stopped the operation (class 7xxx)."""


class MolframeIOError(ResourceError, OSError):
    """An input could not be read or an output could not be written."""


class MemoryBudgetError(ResourceError, MemoryError):
    """The operation needs more memory than the execution budget allows."""


class Cancelled(ResourceError):  # noqa: N818
    """The operation was cancelled before it finished."""


class InternalError(MolframeError, RuntimeError):
    """An internal invariant was violated; this is a bug in MolFrame."""


class MolframeWarning(UserWarning):
    """A non-fatal finding about data."""

    code: str | None
    remedy: str | None

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        remedy: str | None = None,
    ) -> None:
        """Record the message, the diagnostic code and its remedy."""
        super().__init__(message)
        self.code = code
        self.remedy = remedy


class QueryWarning(MolframeWarning):
    """A selection query that is valid but probably not what was meant."""
