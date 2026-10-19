"""Atom-level helpers composed from native kernels.

Each helper is one native call plus array bookkeeping, so it costs what the kernel it
wraps costs. Nothing here loops over atoms in Python, and nothing copies coordinates
beyond the single array the kernel already takes.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._native import Selection, Structure
from .spatial import knn

_NO_ATOM = np.uint32(np.iinfo(np.uint32).max)


@dataclass(frozen=True, slots=True)
class Nearest:
    """The closest atom to each query point.

    ``index`` is the atom's position in the structure; where ``found`` is false it holds
    the sentinel ``4294967295`` and ``distance`` is ``inf``.
    """

    index: NDArray[np.uint32]
    distance: NDArray[np.float32]
    found: NDArray[np.bool_]


def nearest(
    structure: Structure,
    points: ArrayLike,
    *,
    among: Selection | None = None,
    max_distance: float | None = None,
) -> Nearest:
    """Find, for each point, the closest atom of ``structure`` (or of ``among``).

    ``among`` restricts the search to one selection, for example heavy atoms only;
    the returned indices still refer to the whole structure. ``max_distance`` leaves a
    point unfound when no atom lies within that many ångström.
    """
    query = np.ascontiguousarray(points, dtype=np.float32).reshape(-1, 3)
    if among is None:
        atoms: NDArray[np.uint32] | None = None
        coordinates = structure.coordinates
    else:
        if among.is_stale_for(structure):
            message = "selection belongs to a different structure version"
            raise ValueError(message)
        atoms = among.indices
        coordinates = among.to_coordinates()
    if len(coordinates) == 0:
        message = "no atoms to search: the structure or selection is empty"
        raise ValueError(message)
    neighbours, distances = knn(coordinates, 1, query=query, max_distance=max_distance)
    row = neighbours[:, 0]
    found = row != _NO_ATOM
    if atoms is None:
        index = np.where(found, row, _NO_ATOM)
    else:
        index = np.where(found, atoms[np.where(found, row, 0)], _NO_ATOM)
    return Nearest(
        index=index.astype(np.uint32, copy=False),
        distance=distances[:, 0].astype(np.float32, copy=False),
        found=found,
    )


def extract(structure: Structure, selection: Selection) -> Structure:
    """Return a new structure holding only the atoms of ``selection``.

    Bonds whose two ends both survive are kept, as the editor keeps them for any
    deletion; a bond that crosses the boundary is dropped with its removed end.
    """
    if selection.is_stale_for(structure):
        message = "selection belongs to a different structure version"
        raise ValueError(message)
    keep = np.zeros(structure.atom_count, dtype=np.bool_)
    keep[selection.indices] = True
    editor = structure.edit()
    editor.delete(Selection.from_mask(structure, ~keep))
    return editor.finish()


__all__ = ["Nearest", "extract", "nearest"]
