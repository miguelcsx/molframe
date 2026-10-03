//! Typed saccharide metadata, ring geometry and connectivity provenance.

use molframe_core::contract::DictionaryVersion;
use molframe_core::{AtomIndex, BondProvenance, ResidueIndex};

/// SNFG symbols and the generic ring prisms used by molecular renderers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SnfgShape {
    /// Hexose.
    FilledSphere,
    /// N-acetylhexosamine.
    FilledCube,
    /// Hexosamine; crossed secondary colour.
    CrossedCube,
    /// Hexuronate; divided secondary colour.
    DividedDiamond,
    /// Deoxyhexose.
    FilledCone,
    /// N-acetyl-deoxyhexosamine; divided secondary colour.
    DividedCone,
    /// Dideoxyhexose.
    FlatBox,
    /// Pentose.
    FilledStar,
    /// Deoxynonulosonate.
    FilledDiamond,
    /// Dideoxynonulosonate.
    FlatDiamond,
    /// Other or unknown monosaccharide.
    FlatHexagon,
    /// Assigned ketose or apiose symbol.
    Pentagon,
    /// Generic four-membered ring prism.
    DiamondPrism,
    /// Generic five-membered ring prism.
    PentagonalPrism,
    /// Generic six-membered ring prism.
    HexagonalPrism,
    /// Generic seven-membered ring prism.
    HeptagonalPrism,
}

impl SnfgShape {
    /// Stable snake-case name used by the Python surface.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::FilledSphere => "filled_sphere",
            Self::FilledCube => "filled_cube",
            Self::CrossedCube => "crossed_cube",
            Self::DividedDiamond => "divided_diamond",
            Self::FilledCone => "filled_cone",
            Self::DividedCone => "divided_cone",
            Self::FlatBox => "flat_box",
            Self::FilledStar => "filled_star",
            Self::FilledDiamond => "filled_diamond",
            Self::FlatDiamond => "flat_diamond",
            Self::FlatHexagon => "flat_hexagon",
            Self::Pentagon => "pentagon",
            Self::DiamondPrism => "diamond_prism",
            Self::PentagonalPrism => "pentagonal_prism",
            Self::HexagonalPrism => "hexagonal_prism",
            Self::HeptagonalPrism => "heptagonal_prism",
        }
    }
}

/// Curated SNFG identity, independent of a particular deposited residue.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SnfgSymbol {
    /// Community abbreviation, not the CCD identifier.
    pub abbreviation: &'static str,
    /// Human-readable monosaccharide name.
    pub name: &'static str,
    /// Packed sRGB colour, 0xRRGGBB.
    pub color: u32,
    /// Symbol geometry.
    pub shape: SnfgShape,
    /// Packed secondary colour for divided symbols, otherwise absent.
    pub secondary_color: Option<u32>,
}

/// Geometry of a complete finite, nondegenerate ring in ångström.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RingGeometry {
    /// Arithmetic centroid of ring atoms.
    pub center: [f32; 3],
    /// Unit polygon-area normal; sign follows the reported ring atom order.
    pub normal: [f32; 3],
    /// Unit vector from the centroid toward the anomeric carbon, if identified.
    pub anomeric_direction: Option<[f32; 3]>,
}

/// One observed saccharide ring, with no fabricated atoms or geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct Monosaccharide {
    /// Structure-local residue ordinal.
    pub residue: ResidueIndex,
    /// Ring atoms in bonded cyclic order, starting with the oxygen.
    pub ring_atoms: Vec<AtomIndex>,
    /// Topologically identified anomeric carbon; absent if ambiguous.
    pub anomeric_atom: Option<AtomIndex>,
    /// Curated identity, or the generic unknown symbol for an unmapped sugar.
    pub symbol: SnfgSymbol,
    /// Absent for missing/nonfinite coordinates or degenerate geometry.
    pub geometry: Option<RingGeometry>,
}

/// A glycosidic bond directed from its anomeric carbon to the acceptor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CarbohydrateLink {
    /// Donor monosaccharide ordinal in the report.
    pub donor: usize,
    /// Acceptor monosaccharide ordinal, or absent for a terminal attachment.
    pub acceptor: Option<usize>,
    /// Anomeric carbon endpoint.
    pub donor_atom: AtomIndex,
    /// Oxygen, nitrogen or sulfur acceptor endpoint.
    pub acceptor_atom: AtomIndex,
    /// Original bond provenance, or distance inference for the fallback.
    pub provenance: BondProvenance,
}

/// Whether to supplement existing bonds with conservative linkage inference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CarbohydrateOptions {
    /// Search vacant chemically eligible sites within 2 Å; never infer ring edges.
    pub spatial_fallback: bool,
}

impl Default for CarbohydrateOptions {
    fn default() -> Self {
        Self {
            spatial_fallback: true,
        }
    }
}

/// Stable ring and linkage inventory; the input structure is not mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct CarbohydrateReport {
    /// Observed rings in residue and canonical atom order.
    pub monosaccharides: Vec<Monosaccharide>,
    /// Sugar-to-sugar bonds, directed from donor to acceptor.
    pub links: Vec<CarbohydrateLink>,
    /// Sugar-to-nonsugar attachments.
    pub terminal_links: Vec<CarbohydrateLink>,
    /// Recognized sugar residues with no complete observed five/six-member ring.
    pub incomplete_residues: Vec<ResidueIndex>,
    /// Provider version if dictionary topology was requested.
    pub dictionary_version: Option<DictionaryVersion>,
}
