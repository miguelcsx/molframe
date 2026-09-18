//! Public IMD values and explicit resource policy.

use std::time::Duration;

/// Default maximum coordinate or force entries accepted in one message.
pub const DEFAULT_MAX_IMD_ATOMS: usize = 50_000_000;

/// Hard bounds applied before allocating for a peer-provided message.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ImdLimits {
    /// Maximum atoms in one coordinate or force message.
    pub max_atoms: usize,
}

impl Default for ImdLimits {
    fn default() -> Self {
        Self {
            max_atoms: DEFAULT_MAX_IMD_ATOMS,
        }
    }
}

/// TCP connection behavior and protocol resource limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ImdConnectionOptions {
    /// Maximum peer-provided allocation size.
    pub limits: ImdLimits,
    /// Optional socket read timeout.
    pub read_timeout: Option<Duration>,
    /// Optional socket write timeout.
    pub write_timeout: Option<Duration>,
    /// Whether to disable Nagle buffering for latency-sensitive steering.
    pub no_delay: bool,
}

impl Default for ImdConnectionOptions {
    fn default() -> Self {
        Self {
            limits: ImdLimits::default(),
            read_timeout: None,
            write_timeout: None,
            no_delay: true,
        }
    }
}

/// Byte order used by the simulation peer for payload scalars.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImdPeerEndian {
    /// Peer payloads use little-endian scalars.
    Little,
    /// Peer payloads use big-endian scalars.
    Big,
}

/// One complete IMD energy block in canonical molframe units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImdEnergies {
    /// Simulation step reported by the peer.
    pub step: i32,
    /// Temperature in kelvin.
    pub temperature: f32,
    /// Total energy in kJ/mol.
    pub total: f32,
    /// Potential energy in kJ/mol.
    pub potential: f32,
    /// van der Waals energy in kJ/mol.
    pub van_der_waals: f32,
    /// Electrostatic energy in kJ/mol.
    pub electrostatic: f32,
    /// Bond energy in kJ/mol.
    pub bond: f32,
    /// Angle energy in kJ/mol.
    pub angle: f32,
    /// Dihedral energy in kJ/mol.
    pub dihedral: f32,
    /// Improper energy in kJ/mol.
    pub improper: f32,
}

/// One atom-targeted steering force in canonical kJ mol⁻¹ Å⁻¹.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImdForce {
    /// Zero-based simulation atom index.
    pub atom: u32,
    /// Cartesian force vector.
    pub force: [f32; 3],
}

/// One complete message received from an IMD peer.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum ImdMessage {
    /// The simulation requested a clean detach while continuing to run.
    Disconnect,
    /// A complete energy block.
    Energies(ImdEnergies),
    /// Coordinates in ångström, in simulation atom order.
    Coordinates(Vec<[f32; 3]>),
    /// Start or resume control message.
    Go,
    /// The simulation requested termination of the session and job.
    Kill,
    /// Atom-targeted forces, converted to canonical units.
    Forces(Vec<ImdForce>),
    /// Pause-toggle control message.
    Pause,
    /// Updated coordinate transmission interval in simulation steps.
    TransmissionRate(u32),
}
