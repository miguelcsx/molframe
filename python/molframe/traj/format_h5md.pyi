from .._trajectory import Timestep

class H5mdError(Exception): ...
class H5mdUnitSystem:
    def __init__(self, length_unit: str, length_to_angstrom: float, time_unit: str, time_to_picosecond: float, velocity_unit: str, velocity_to_angstrom_per_picosecond: float, force_unit: str, force_to_kilojoule_per_mole_angstrom: float) -> None: ...
    @staticmethod
    def canonical() -> H5mdUnitSystem: ...
    length_unit: str
    length_to_angstrom: float
    time_unit: str
    time_to_picosecond: float
    velocity_unit: str
    velocity_to_angstrom_per_picosecond: float
    force_unit: str
    force_to_kilojoule_per_mole_angstrom: float
class H5mdOptions:
    def __init__(self, *, particle_group: str = ..., units: H5mdUnitSystem | None = ...) -> None: ...
    particle_group: str
    units: H5mdUnitSystem
class H5mdMetadata:
    options: H5mdOptions
    creator_name: str
    creator_version: str
class H5mdTrajectory:
    frames: list[Timestep]
    metadata: H5mdMetadata
def parse_h5md(bytes: bytes) -> list[Timestep]: ...
def parse_h5md_with_options(bytes: bytes, options: H5mdOptions) -> list[Timestep]: ...
def parse_h5md_record_with_options(bytes: bytes, options: H5mdOptions) -> H5mdTrajectory: ...
def write_h5md(frames: list[Timestep], options: H5mdOptions) -> bytes: ...
def write_h5md_with_metadata(frames: list[Timestep], metadata: H5mdMetadata) -> bytes: ...
