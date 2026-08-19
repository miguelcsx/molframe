"""Trajectory format records, streaming controls, and intrinsic value types."""

from os import PathLike
from typing import final
from numpy import float32, float64, int64, uint32
from numpy.typing import NDArray
from numpy.ma import MaskedArray
from .surface import SurfaceGridOptions
from ._native import CartesianFit
from ._trajectory_dispatch import AmberAsciiReadOptions, TrajectoryReadOptions, read_trajectory, write_trajectory
DEFAULT_PAIRWISE_MEMORY_LIMIT: int

@final
class ImdLimits:
    def __init__(self, max_atoms: int | None = ...) -> None: ...
    max_atoms: int

@final
class ImdConnectionOptions:
    def __init__(self, limits: ImdLimits | None = ..., *, read_timeout: float | None = ..., write_timeout: float | None = ..., no_delay: bool = ...) -> None: ...
    @property
    def limits(self) -> ImdLimits: ...
    @property
    def read_timeout(self) -> float | None: ...
    @property
    def write_timeout(self) -> float | None: ...
    @property
    def no_delay(self) -> bool: ...

@final
class ImdPeerEndian:
    Little: ImdPeerEndian
    Big: ImdPeerEndian

@final
class ImdEnergies:
    def __init__(self, step: int, temperature: float, total: float, potential: float, van_der_waals: float, electrostatic: float, bond: float, angle: float, dihedral: float, improper: float) -> None: ...
    step: int
    temperature: float
    total: float
    potential: float
    van_der_waals: float
    electrostatic: float
    bond: float
    angle: float
    dihedral: float
    improper: float

@final
class ImdForce:
    def __init__(self, atom: int, force: tuple[float, float, float]) -> None: ...
    atom: int
    force: tuple[float, float, float]

@final
class ImdMessage:
    @staticmethod
    def from_disconnect() -> ImdMessage: ...
    @staticmethod
    def from_energies(value: ImdEnergies) -> ImdMessage: ...
    @staticmethod
    def from_coordinates(value: list[tuple[float, float, float]]) -> ImdMessage: ...
    @staticmethod
    def from_go() -> ImdMessage: ...
    @staticmethod
    def from_kill() -> ImdMessage: ...
    @staticmethod
    def from_forces(value: list[ImdForce]) -> ImdMessage: ...
    @staticmethod
    def from_pause() -> ImdMessage: ...
    @staticmethod
    def from_transmission_rate(rate: int) -> ImdMessage: ...
    def __repr__(self) -> str: ...
    kind: str
    energies: ImdEnergies | None
    coordinates: list[tuple[float, float, float]] | None
    forces: list[ImdForce] | None
    rate: int | None

@final
class ImdClient:
    @staticmethod
    def connect(address: str, *, options: ImdConnectionOptions | None = ...) -> ImdClient: ...
    def receive(self) -> ImdMessage: ...
    def send_forces(self, forces: list[ImdForce]) -> None: ...
    def toggle_pause(self) -> None: ...
    def resume(self) -> None: ...
    def set_transmission_rate(self, steps: int) -> None: ...
    def disconnect(self) -> None: ...
    def kill(self) -> None: ...
    @property
    def peer_endian(self) -> ImdPeerEndian: ...
    @property
    def closed(self) -> bool: ...
from ._io_types import DielectricError, DmsError, EnsembleGeometryError, EnsembleSimilarityError, EnsembleStatisticsError, GovernedEnsembleError, ImdError, KMeansError, MsdError, PathSimilarityError, TrajectoryError, TrajectoryIoError, UpdatingSelectionError, WaterDynamicsError
from ._trajectory_dms import DmsBond, DmsCell, DmsFrame, DmsParticle, DmsSystem, DmsTopology, DmsVersion, read_dms, write_dms
from ._trajectory_kernels import DielectricEstimate, MeanSquaredDisplacement, TrajectoryDielectricOptions, dielectric_from_dipoles, generalized_procrustes_mean, kmeans, mean_squared_displacement, pairwise_fitted_rmsd, pairwise_torus_distance, rmsd_to_reference
from ._trajectory_ensemble import (
    Clustering, ConvergenceBlock, EnsembleDistanceMatrix, FrameAlignment, GroupVariance,
    HarmonicSimilarity, HarmonicSimilarityOptions, KMeans, KMeansOptions, Linkage,
    RemainderPolicy, analyse_agglomerative_clustering, analyse_block_convergence,
    analyse_cluster_population_similarity, analyse_generalized_procrustes_mean,
    analyse_group_coordinate_variance, analyse_harmonic_ensemble_similarity, analyse_kmeans,
    analyse_medoid, analyse_pairwise_fitted_rmsd, analyse_pairwise_torus_distance,
    analyse_rmsd_to_reference,
)
from .traj.format_amber_netcdf import AmberNetcdfError, AmberNetcdfMetadata, AmberNetcdfPrecision, AmberNetcdfTrajectory, AmberNetcdfWriteOptions, parse_amber_netcdf, parse_amber_netcdf_record, write_amber_netcdf
from .traj.format_gromacs_itp import GromacsInteraction, GromacsItp, GromacsItpAtom, GromacsItpError, GromacsMoleculeType, parse_gromacs_itp
from .traj.format_amber_topology import AmberSection, AmberTopology, AmberTopologyAngle, AmberTopologyAtom, AmberTopologyBond, AmberTopologyDihedral, AmberTopologyError, AmberTopologyResidue, parse_amber_topology
from .traj.format_psf import PsfAtom, PsfError, PsfTopology, parse_psf, write_psf
from .traj.format_namd import NamdBinary, NamdEndian, NamdError, parse_namd_binary, write_namd_binary
from .traj.format_dcd import DcdEndian, DcdError, DcdHeader, DcdTrajectory, DcdWriteOptions, parse_dcd, write_dcd
from .traj.format_trr import TrrError, TrrPrecision, TrrTrajectory, TrrWriteOptions, parse_trr, write_trr, write_trr_with_precisions
from .traj.format_xtc import XtcError, XtcTrajectory, XtcWriteOptions, parse_xtc, write_xtc, write_xtc_with_precisions
from .traj.format_gsd import GsdError, GsdOptions, GsdTrajectory, parse_gsd, write_gsd
from .traj.format_h5md import H5mdError, H5mdMetadata, H5mdOptions, H5mdTrajectory, H5mdUnitSystem, parse_h5md, parse_h5md_record_with_options, parse_h5md_with_options, write_h5md, write_h5md_with_metadata
from .traj.format_tng import TngCompression, TngError, TngTrajectory, TngWriteOptions, parse_tng, write_tng
from .traj.format_trz import TrzError, TrzTrajectory, parse_trz, write_trz
from .traj.format_gromos import GromosBoundary, GromosError, GromosTrajectory, parse_gromos11_trc
from .traj.format_hoomd_xml import HoomdBox, HoomdConfiguration, HoomdInteraction, HoomdXmlError, parse_hoomd_xml
from .traj.format_lammps import LammpsError, parse_lammps_dump
from .traj.format_lammps_data import LammpsAtomStyle, LammpsData, LammpsDataAtom, LammpsDataCell, LammpsDataError, LammpsInteraction, parse_lammps_data
from .traj.format_tpr import TprAtom, TprBond, TprError, TprHeader, TprResidue, TprTopology, parse_tpr

@final
class AimsAtom:
    def __init__(self, species: str, position: tuple[float, float, float]) -> None: ...
    species: str
    position: tuple[float, float, float]
@final
class AimsGeometry:
    def __init__(self, atoms: list[AimsAtom], lattice_vectors: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None = ...) -> None: ...
    atoms: list[AimsAtom]
    lattice_vectors: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None
class AimsError(Exception): ...
def parse_aims_geometry(text: str) -> AimsGeometry: ...
def write_aims_geometry(geometry: AimsGeometry) -> str: ...

@final
class GroAtom:
    def __init__(self, residue_number: int, residue_name: str, atom_name: str, atom_number: int, position: tuple[float, float, float], velocity: tuple[float, float, float] | None = ...) -> None: ...
    residue_number: int
    residue_name: str
    atom_name: str
    atom_number: int
    position: tuple[float, float, float]
    velocity: tuple[float, float, float] | None
@final
class GroFrame:
    def __init__(self, title: str, atoms: list[GroAtom], box_values: list[float]) -> None: ...
    title: str
    atoms: list[GroAtom]
    box_values: list[float]
class GroError(Exception): ...
def parse_gro_records(text: str) -> list[GroFrame]: ...
def write_gro(frames: list[GroFrame]) -> str: ...

@final
class TxyzAtom:
    def __init__(self, id: int, name: str, position: tuple[float, float, float], atom_type: str, bonds: list[int]) -> None: ...
    id: int
    name: str
    position: tuple[float, float, float]
    atom_type: str
    bonds: list[int]
@final
class TxyzFrame:
    def __init__(self, title: str, atoms: list[TxyzAtom]) -> None: ...
    title: str
    atoms: list[TxyzAtom]
class TxyzError(Exception): ...
def parse_txyz_records(text: str) -> list[TxyzFrame]: ...
def write_txyz(frames: list[TxyzFrame]) -> str: ...

@final
class CharmmCardFormat:
    Standard: CharmmCardFormat
    Extended: CharmmCardFormat
    Free: CharmmCardFormat
@final
class CharmmAtom:
    def __init__(self, atom_number: int, residue_number: int, residue_name: str, atom_name: str, position: tuple[float, float, float], segment_id: str, residue_id: str, weight: float) -> None: ...
    atom_number: int
    residue_number: int
    residue_name: str
    atom_name: str
    position: tuple[float, float, float]
    segment_id: str
    residue_id: str
    weight: float
@final
class CharmmCard:
    def __init__(self, titles: list[str], format: CharmmCardFormat, atoms: list[CharmmAtom]) -> None: ...
    titles: list[str]
    format: CharmmCardFormat
    atoms: list[CharmmAtom]
class CharmmError(Exception): ...
def parse_charmm_record(text: str) -> CharmmCard: ...
def write_charmm_card(card: CharmmCard) -> str: ...

@final
class GamessRunType:
    Optimize: GamessRunType
    Surface: GamessRunType
@final
class GamessAtom:
    def __init__(self, name: str, nuclear_charge: float | None = ...) -> None: ...
    name: str
    nuclear_charge: float | None
@final
class GamessFrame:
    def __init__(self, step: int, positions: list[tuple[float, float, float]], energy: float | None = ..., surface_coordinates: tuple[float, float] | None = ...) -> None: ...
    step: int
    energy: float | None
    surface_coordinates: tuple[float, float] | None
    positions: list[tuple[float, float, float]]
@final
class GamessTrajectory:
    def __init__(self, run_type: GamessRunType, atoms: list[GamessAtom], frames: list[GamessFrame]) -> None: ...
    run_type: GamessRunType
    atoms: list[GamessAtom]
    frames: list[GamessFrame]
class GamessError(Exception): ...
def parse_gamess_output(text: str) -> GamessTrajectory: ...

@final
class FrameValue:
    kind: str
    float_value: float | None
    integer_value: int | None
    text_value: str | None
    floats_value: list[float] | None
    @staticmethod
    def float(value: float) -> FrameValue: ...
    @staticmethod
    def integer(value: int) -> FrameValue: ...
    @staticmethod
    def text(value: str) -> FrameValue: ...
    @staticmethod
    def floats(value: list[float]) -> FrameValue: ...

@final
class Timestep:
    def __init__(
        self,
        positions: list[tuple[float, float, float]],
        *,
        frame: int = 0,
        time: float | None = None,
        dt: float | None = None,
        velocities: list[tuple[float, float, float]] | None = None,
        forces: list[tuple[float, float, float]] | None = None,
        cell: object | None = None,
        data: list[tuple[str, FrameValue]] | None = None,
    ) -> None: ...
    frame: int
    time: float | None
    dt: float | None
    positions: list[tuple[float, float, float]]
    velocities: list[tuple[float, float, float]] | None
    forces: list[tuple[float, float, float]] | None
    cell: object | None
    data: list[tuple[str, FrameValue]]

@final
class AmberRestartLayout:
    Auto: AmberRestartLayout
    Coordinates: AmberRestartLayout
    CoordinatesBox3: AmberRestartLayout
    CoordinatesBox6: AmberRestartLayout
    CoordinatesVelocities: AmberRestartLayout
    CoordinatesVelocitiesBox3: AmberRestartLayout
    CoordinatesVelocitiesBox6: AmberRestartLayout

@final
class AmberRestart:
    def __init__(self, title: str, layout: AmberRestartLayout, timestep: Timestep) -> None: ...
    title: str
    layout: AmberRestartLayout
    timestep: Timestep

class AmberError(Exception): ...
def parse_amber_restart(text: str, layout: AmberRestartLayout) -> Timestep: ...
def parse_amber_restart_record(text: str, layout: AmberRestartLayout) -> AmberRestart: ...
def write_amber_restart(record: AmberRestart) -> str: ...
def parse_amber_ascii_trajectory(text: str, atom_count: int, periodic_box: bool) -> list[Timestep]: ...

@final
class DlPolyAtom:
    def __init__(self, name: str, index: int | None = ..., mass: float | None = ..., charge: float | None = ...) -> None: ...
    name: str
    index: int | None
    mass: float | None
    charge: float | None

@final
class DlPolyFrame:
    def __init__(self, step: int, time: float, positions: list[tuple[float, float, float]], velocities: list[tuple[float, float, float]] | None = ..., forces: list[tuple[float, float, float]] | None = ..., lattice_vectors: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None = ...) -> None: ...
    step: int
    time: float
    positions: list[tuple[float, float, float]]
    velocities: list[tuple[float, float, float]] | None
    forces: list[tuple[float, float, float]] | None
    lattice_vectors: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None
    def to_timestep(self, frame: int) -> Timestep: ...

@final
class DlPolyConfig:
    def __init__(self, title: str, level: int, boundary: int, atoms: list[DlPolyAtom], frame: DlPolyFrame) -> None: ...
    title: str
    level: int
    boundary: int
    atoms: list[DlPolyAtom]
    frame: DlPolyFrame

@final
class DlPolyHistory:
    def __init__(self, title: str, level: int, boundary: int, atoms: list[DlPolyAtom], frames: list[DlPolyFrame]) -> None: ...
    title: str
    level: int
    boundary: int
    atoms: list[DlPolyAtom]
    frames: list[DlPolyFrame]

class DlPolyError(Exception): ...
def parse_dlpoly_config(text: str) -> DlPolyConfig: ...
def parse_dlpoly_history(text: str) -> DlPolyHistory: ...
def write_dlpoly_config(config: DlPolyConfig) -> str: ...
def write_dlpoly_history(history: DlPolyHistory) -> str: ...

@final
class XyzAtom:
    def __init__(self, element: str, position: tuple[float, float, float]) -> None: ...
    element: str
    position: tuple[float, float, float]

@final
class XyzFrame:
    def __init__(self, comment: str, atoms: list[XyzAtom]) -> None: ...
    comment: str
    atoms: list[XyzAtom]

def parse_xyz(text: str) -> list[XyzFrame] | None: ...
def write_xyz(frames: list[XyzFrame]) -> str: ...


@final
class PeriodicAngle:
    def __init__(self, radians: float) -> None: ...
    @staticmethod
    def from_radians(radians: float) -> PeriodicAngle: ...
    @staticmethod
    def from_degrees(degrees: float) -> PeriodicAngle: ...
    @property
    def radians(self) -> float: ...
    @property
    def degrees(self) -> float: ...
    def delta_to(self, other: PeriodicAngle) -> float: ...
    def signed_delta(self, other: PeriodicAngle) -> float: ...
    def distance_to(self, other: PeriodicAngle) -> float: ...
    def distance(self, other: PeriodicAngle) -> float: ...

@final
class Rotation3:
    IDENTITY: Rotation3
    def __init__(self, matrix: list[list[float]]) -> None: ...
    @staticmethod
    def identity() -> Rotation3: ...
    @staticmethod
    def from_matrix(matrix: list[list[float]], tolerance: float) -> Rotation3: ...
    @staticmethod
    def exponential(rotation_vector: list[float]) -> Rotation3: ...
    @staticmethod
    def exp(rotation_vector: list[float]) -> Rotation3: ...
    @staticmethod
    def exponential_with_options(rotation_vector: list[float], options: object) -> Rotation3: ...
    @staticmethod
    def exp_with_options(rotation_vector: list[float], options: object) -> Rotation3: ...
    @property
    def matrix(self) -> list[list[float]]: ...
    def logarithm(self) -> list[float]: ...
    def log(self) -> list[float]: ...
    def log_with_options(self, options: object) -> list[float]: ...
    def apply(self, vector: list[float]) -> list[float]: ...
    def inverse(self) -> Rotation3: ...
    def then(self, next: Rotation3) -> Rotation3: ...
    def distance_to(self, other: Rotation3) -> float: ...
    def interpolate(self, other: Rotation3, fraction: float) -> Rotation3: ...
    def interpolate_with_options(self, other: Rotation3, fraction: float, options: object) -> Rotation3: ...

@final
class SurfaceMesh:
    @property
    def vertices(self) -> NDArray[float32]: ...
    @property
    def faces(self) -> NDArray[uint32]: ...
    def curvature(self) -> list[tuple[float, float, float, float, float, float, str]]: ...

@final
class PcaResult:
    mean: NDArray[float64]
    eigenvalues: NDArray[float64]
    components: NDArray[float64]
    projections: NDArray[float64]

@final
class DiffusionMap:
    eigenvalues: NDArray[float64]
    coordinates: NDArray[float64]
    graph_components: NDArray[int64]

def surface_mesh(positions: NDArray[float32], radii: NDArray[float32], probe: float, options: SurfaceGridOptions) -> SurfaceMesh: ...
