"""Dividing work so the answer does not depend on how many threads ran it."""

from .._native import (
    DEFAULT_BLOCK_ITEMS,
    BlockPlan,
    try_for_each_block_in,
    ReductionPolicy,
)

__all__ = [
    "DEFAULT_BLOCK_ITEMS",
    "BlockPlan",
    "try_for_each_block_in",
    "ReductionPolicy",
]
