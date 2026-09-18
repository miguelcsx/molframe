from typing import final
from .._trajectory import Timestep

@final
class NamdEndian:
    Little: NamdEndian
    Big: NamdEndian

@final
class NamdBinary:
    frame: Timestep
    endian: NamdEndian

class NamdError(Exception): ...
def parse_namd_binary(bytes: bytes) -> NamdBinary: ...
def write_namd_binary(frame: Timestep, endian: NamdEndian) -> bytes: ...

__all__: list[str]
