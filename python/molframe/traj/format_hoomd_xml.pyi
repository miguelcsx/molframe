class HoomdXmlError(Exception): ...
class HoomdBox:
    lengths: tuple[float, float, float]
    tilt: tuple[float, float, float]
class HoomdInteraction:
    kind: str
    particles: list[int]
class HoomdConfiguration:
    step: int | None
    dimensions: int | None
    particle_count: int
    cell: HoomdBox | None
    positions: list[tuple[float, float, float]]
    images: list[tuple[int, int, int]]
    velocities: list[tuple[float, float, float]]
    accelerations: list[tuple[float, float, float]]
    orientations: list[tuple[float, float, float, float]]
    types: list[str]
    masses: list[float]
    charges: list[float]
    diameters: list[float]
    bodies: list[int]
    bonds: list[HoomdInteraction]
    angles: list[HoomdInteraction]
    dihedrals: list[HoomdInteraction]
    impropers: list[HoomdInteraction]
    extensions: dict[str, str]
def parse_hoomd_xml(source: str) -> HoomdConfiguration: ...

__all__: list[str]
