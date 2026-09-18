from .._trajectory import Timestep

class GromosBoundary:
    Vacuum: GromosBoundary
    Rectangular: GromosBoundary
    Triclinic: GromosBoundary
    TruncatedOctahedron: GromosBoundary
class GromosError(Exception): ...
class GromosTrajectory:
    title: str
    frames: list[Timestep]
    boundaries: list[GromosBoundary | None]
def parse_gromos11_trc(source: str) -> GromosTrajectory: ...
