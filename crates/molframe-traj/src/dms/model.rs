//! Lossless typed values for the standard DMS chemical-structure tables.

/// DMS schema version when the source declares one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DmsVersion {
    /// Major schema version.
    pub major: u32,
    /// Minor schema version.
    pub minor: u32,
}

/// One particle's identity and standard chemical annotations.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DmsParticle {
    /// Atomic number; zero remains meaningful for pseudo-particles.
    pub atomic_number: Option<u16>,
    /// Component identifier.
    pub component: Option<i64>,
    /// Non-bonded parameter type.
    pub nonbonded_type: Option<i64>,
    /// Mass in atomic mass units.
    pub mass: Option<f64>,
    /// Partial charge in elementary-charge units.
    pub charge: Option<f64>,
    /// Residue number.
    pub residue_id: Option<i64>,
    /// Residue name.
    pub residue_name: Option<Box<str>>,
    /// Chain identifier.
    pub chain: Option<Box<str>>,
    /// Segment identifier.
    pub segment: Option<Box<str>>,
    /// Atom or pseudo-particle name.
    pub name: Option<Box<str>>,
    /// Residue insertion code.
    pub insertion: Option<Box<str>>,
    /// Formal charge in elementary-charge units.
    pub formal_charge: Option<f64>,
    /// PDB occupancy.
    pub occupancy: Option<f64>,
    /// PDB temperature factor.
    pub b_factor: Option<f64>,
    /// Temperature group identifier.
    pub temperature_group: Option<i64>,
    /// Energy group identifier.
    pub energy_group: Option<i64>,
    /// Ligand group identifier.
    pub ligand_group: Option<i64>,
    /// Force-bias group identifier.
    pub bias_group: Option<i64>,
}

/// A chemical bond between zero-based particle identifiers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmsBond {
    /// Lower particle identifier.
    pub p0: u32,
    /// Higher particle identifier.
    pub p1: u32,
    /// Bond order as stored by DMS.
    pub order: f64,
}

/// Triclinic periodic cell vectors in Angstrom.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmsCell {
    /// Cell vectors in file order.
    pub vectors: [[f64; 3]; 3],
}

/// The single coordinate frame carried by a DMS structure database.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DmsFrame {
    /// Positions in Angstrom, ordered by particle identifier.
    pub positions: Vec<[f64; 3]>,
    /// Optional velocities in Angstrom per picosecond, ordered by particle id.
    pub velocities: Vec<Option<[f64; 3]>>,
    /// Optional periodic cell.
    pub cell: Option<DmsCell>,
}

/// Chemical topology carried by the DMS structural tables.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DmsTopology {
    /// Particles ordered by their implicit, contiguous identifier.
    pub particles: Vec<DmsParticle>,
    /// Bonds in ascending `(p0, p1)` order.
    pub bonds: Vec<DmsBond>,
}

/// Complete standard DMS chemical structure and its coordinate frame.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DmsSystem {
    /// Optional declared schema version.
    pub version: Option<DmsVersion>,
    /// Particle and bond topology.
    pub topology: DmsTopology,
    /// Coordinates, optional velocities, and optional cell.
    pub frame: DmsFrame,
}
