"""The compiled extension module, ``molframe._native``.

Every name the extension binds is declared once, by the stub module that owns
it, so this gathers the two halves of the core vocabulary and adds what only the
flat module carries. The namespaces the extension exposes as submodules have
their own stub beside this one.
"""

from .._core_primitives import *
from .._core_topology import *
from .validate import classify_ramachandran as classify_ramachandran
