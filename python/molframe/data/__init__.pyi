import builtins
from os import PathLike
from typing import final

@final
class DataEntry:
    """One manifest entry and the state of its cache file."""

    @property
    def name(self) -> str: ...
    @property
    def version(self) -> str: ...
    @property
    def licence(self) -> str: ...
    @property
    def path(self) -> str: ...
    @property
    def status(self) -> str:
        """``missing``, ``current``, ``stale`` or ``invalid``."""
    @property
    def detail(self) -> str | None:
        """Why the entry is stale or invalid; ``None`` otherwise."""

@final
class DownloadLimits:
    """The explicit resource bounds of one download."""

    def __init__(self, *, max_bytes: int, timeout_seconds: int, redirect_limit: int) -> None: ...
    @property
    def max_bytes(self) -> int: ...
    @property
    def timeout_seconds(self) -> int: ...
    @property
    def redirect_limit(self) -> int: ...

@final
class Fetched:
    """A table that ``fetch`` just stored."""

    @property
    def path(self) -> str: ...
    @property
    def bytes(self) -> int: ...
    @property
    def rows(self) -> int: ...

def list(
    manifest: str | PathLike[str], directory: str | PathLike[str] | None = None
) -> builtins.list[DataEntry]:
    """Return every manifest entry with the state of its cache file; no network.

    ``directory`` defaults to the table store's cache directory.
    """

def verify(
    manifest: str | PathLike[str],
    directory: str | PathLike[str] | None = None,
    *,
    name: str | None = None,
) -> builtins.list[DataEntry]:
    """Re-read the cache files of the manifest (or the one ``name``); no network.

    A malformed manifest raises ``E2201``; an unknown ``name`` raises ``E7102``.
    """

def fetch(
    manifest: str | PathLike[str],
    name: str,
    limits: DownloadLimits,
    directory: str | PathLike[str] | None = None,
    *,
    replace: bool = False,
) -> Fetched:
    """Download ``name``, convert, validate and store it; the only call that uses the network.

    The bytes are accepted only when their SHA-256 matches the manifest. An existing cache
    file raises ``E7901`` unless ``replace`` is true.
    """
