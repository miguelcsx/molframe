//! Numerical and geometric definition shared by automatic and explicit DSSP.

use std::ops::RangeInclusive;

/// Numerical and pattern parameters for the native DSSP 4 classifier.
#[derive(Clone, Debug, PartialEq)]
pub struct DsspOptions {
    /// Electrostatic prefactor in kcal/mol angstroms.
    pub electrostatic_prefactor: f64,
    /// Interactions below this energy are hydrogen bonds.
    pub hydrogen_bond_energy: f64,
    /// Reconstructed amide N-H distance in angstroms.
    pub amide_hydrogen_distance: f32,
    /// Minimum index separation; 1 includes the standard nearest-neighbour search.
    pub minimum_sequence_separation: usize,
    /// Offset of consecutive alpha-helical turns.
    pub helix_offset: usize,
    /// Offset of consecutive three-ten-helical turns.
    pub three_ten_offset: usize,
    /// Offset of consecutive pi-helical turns.
    pub pi_offset: usize,
    /// Inclusive offsets identifying hydrogen-bonded turns.
    pub turn_offsets: RangeInclusive<usize>,
    /// Minimum change of C-alpha direction for a bend, in degrees.
    pub bend_angle_degrees: f32,
}

impl Default for DsspOptions {
    fn default() -> Self {
        Self {
            electrostatic_prefactor: 27.888,
            hydrogen_bond_energy: -0.5,
            amide_hydrogen_distance: 1.0,
            minimum_sequence_separation: 1,
            helix_offset: 4,
            three_ten_offset: 3,
            pi_offset: 5,
            turn_offsets: 3..=5,
            bend_angle_degrees: 70.0,
        }
    }
}

impl DsspOptions {
    /// Whether every parameter is finite and has a meaningful domain.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.electrostatic_prefactor.is_finite()
            && self.electrostatic_prefactor > 0.0
            && self.hydrogen_bond_energy.is_finite()
            && self.amide_hydrogen_distance.is_finite()
            && self.amide_hydrogen_distance > 0.0
            && self.minimum_sequence_separation > 0
            && self.helix_offset > self.minimum_sequence_separation
            && self.three_ten_offset > self.minimum_sequence_separation
            && self.pi_offset > self.minimum_sequence_separation
            && !self.turn_offsets.is_empty()
            && *self.turn_offsets.start() > self.minimum_sequence_separation
            && self.bend_angle_degrees.is_finite()
            && (0.0..=180.0).contains(&self.bend_angle_degrees)
    }
}

/// A native secondary-structure parameter is nonfinite or outside its domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InvalidDsspOptions;
impl std::fmt::Display for InvalidDsspOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("secondary-structure options are invalid")
    }
}
impl std::error::Error for InvalidDsspOptions {}

/// One residue's backbone, supplied by either atom-name or semantic-role extraction.
#[derive(Clone, Copy, Debug, Default)]
pub struct DsspBackbone {
    /// Model identity; hydrogen bonds never connect different models.
    pub model: u32,
    /// Chain identity; local patterns never cross chain boundaries.
    pub chain: u32,
    /// Proline cannot donate a backbone hydrogen bond.
    pub proline: bool,
    /// Alpha-carbon position in angstroms.
    pub ca: Option<[f32; 3]>,
    /// Carbonyl-carbon position in angstroms.
    pub carbon: Option<[f32; 3]>,
    /// Carbonyl-oxygen position in angstroms.
    pub oxygen: Option<[f32; 3]>,
    /// Backbone-nitrogen position in angstroms.
    pub nitrogen: Option<[f32; 3]>,
}

impl DsspBackbone {
    /// Complete heavy backbone, including oxygen; donation is not required.
    #[must_use]
    pub const fn is_evaluable(&self) -> bool {
        self.ca.is_some()
            && self.nitrogen.is_some()
            && self.carbon.is_some()
            && self.oxygen.is_some()
    }
}
