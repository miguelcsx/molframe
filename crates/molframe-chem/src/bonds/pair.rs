//! Pairwise bonding predicates and per-element reach used by bond perception.
//!
//! Every test here is constant time and works on squared distances, so the
//! grid search stays free of square roots; the thresholds are reference data
//! for which elements may bond inferentially at all.

use molframe_core::structure::AtomRef;

/// Whether a distance lies between the closest plausible contact and the sum of
/// the two elements' bonding thresholds.
pub(super) fn in_bonding_window(
    left_threshold: f32,
    right_threshold: f32,
    distance_squared: f32,
) -> bool {
    let maximum = maximum_bonding_distance(left_threshold, right_threshold);
    distance_squared >= 0.16 && distance_squared <= maximum * maximum
}

pub(super) fn maximum_bonding_distance(left_threshold: f32, right_threshold: f32) -> f32 {
    (left_threshold + right_threshold) / 1.95
}

/// Whether two atoms may bond at all: two hydrogens never do, and two atoms in
/// different alternate locations do not coexist.
pub(super) fn pair_compatible(left: AtomRef<'_>, right: AtomRef<'_>) -> bool {
    let both_hydrogen = left
        .element()
        .is_some_and(molframe_core::Element::is_hydrogen)
        && right
            .element()
            .is_some_and(molframe_core::Element::is_hydrogen);
    !both_hydrogen && alt_compatible(left, right)
}

pub(super) fn alt_compatible(left: AtomRef<'_>, right: AtomRef<'_>) -> bool {
    match (left.alt_id(), right.alt_id()) {
        (Some(left), Some(right)) => left.is_blank() || right.is_blank() || left == right,
        _ => true,
    }
}

pub(super) fn within(left: AtomRef<'_>, right: AtomRef<'_>, maximum: f32) -> bool {
    let (Some(left), Some(right)) = (left.position(), right.position()) else {
        return false;
    };
    squared_distance(left, right) <= maximum * maximum
}

pub(super) fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}

/// The bonding reach of one element, in ångström, or `None` when this reader
/// has no reference value for it.
///
/// An element with no entry takes part in no inferred bond at all: perception
/// leaves it to deposited connectivity rather than guessing a radius. Selenium
/// is the case that matters in practice — selenomethionine's `SE` is a real
/// covalent substituent — but the same list governs every other element, and
/// widening it is a data decision, not a code one.
pub(super) fn element_threshold(element: molframe_core::Element) -> Option<f32> {
    match element.symbol() {
        "H" => Some(1.42),
        "C" => Some(1.75),
        "N" => Some(1.60),
        "O" => Some(1.52),
        "P" => Some(2.00),
        "S" => Some(1.90),
        "F" | "Cl" | "Br" | "I" => Some(1.80),
        symbol if is_metal(symbol) => metal_threshold(symbol),
        _ => None,
    }
}

/// The coordination reach of one metal, in ångström, on the same scale as the
/// non-metal thresholds.
///
/// A metal–ligand bond is accepted up to `(metal + ligand) / 1.95`, so a metal's
/// value is chosen so that the longest ordinary coordination distance of that
/// metal to a nitrogen, oxygen or sulfur donor still falls inside the window. One
/// shared value cannot do it: a heme iron sits 1.9–2.3 Å from its pyrrole
/// nitrogens and a calcium 2.3–2.6 Å from its oxygens, so a single compromise
/// either drops iron–nitrogen bonds that are really there or lets small metals
/// reach through to neighbours that are not coordinated.
fn metal_threshold(symbol: &str) -> Option<f32> {
    match symbol {
        "Li" | "Al" => Some(2.80),
        "Na" => Some(3.30),
        "K" => Some(3.90),
        "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Au" | "Cr" | "V" => Some(2.90),
        "Mg" | "Hg" => Some(3.00),
        "Ca" => Some(3.40),
        "Zn" => Some(2.85),
        "Cd" | "Mo" | "W" | "Ag" => Some(3.10),
        _ => None,
    }
}

pub(super) fn is_metal(symbol: &str) -> bool {
    matches!(
        symbol,
        "Li" | "Na"
            | "K"
            | "Mg"
            | "Ca"
            | "Mn"
            | "Fe"
            | "Co"
            | "Ni"
            | "Cu"
            | "Zn"
            | "Mo"
            | "Ag"
            | "Cd"
            | "Hg"
            | "Au"
            | "Al"
            | "Cr"
            | "V"
            | "W"
    )
}
