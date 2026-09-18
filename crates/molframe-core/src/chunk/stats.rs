//! What a chunk can say about itself without being read.
//!
//! Every chunk carries a small summary, and the summary is the single largest
//! source of speed in the whole design. Asked for zinc atoms within five
//! ångström of one chain, the engine reads no rows at all from a chunk whose
//! element set lacks zinc, and none from a chunk whose bounding box cannot come
//! that close. On an ordinary entry that removes well over nine rows in ten
//! before the first comparison.
//!
//! The summary is built once, when the chunk is, and is a few dozen bytes
//! against a quarter of a megabyte of atom data.

use crate::coords::Aabb;
use crate::element::Element;

/// The set of elements present, one bit per atomic number.
///
/// The periodic table fits in two machine words, so membership is a shift and a
/// test rather than a search.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
#[repr(transparent)]
pub struct ElementMask(u128);

impl ElementMask {
    /// The empty set.
    pub const EMPTY: Self = Self(0);

    /// Adds an element.
    pub const fn insert(&mut self, element: Element) {
        self.0 |= 1u128 << element.atomic_number();
    }

    /// Returns true when the element may be present.
    #[must_use]
    pub const fn contains(self, element: Element) -> bool {
        self.0 >> element.atomic_number() & 1 == 1
    }

    /// Returns true when the set is empty.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// Adds everything in `other`.
    pub const fn union_with(&mut self, other: Self) {
        self.0 |= other.0;
    }

    /// Returns true when the two sets share an element.
    #[must_use]
    pub const fn intersects(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }

    /// The number of distinct elements.
    #[must_use]
    pub const fn len(self) -> u32 {
        self.0.count_ones()
    }
}

impl FromIterator<Element> for ElementMask {
    fn from_iter<T: IntoIterator<Item = Element>>(iter: T) -> Self {
        let mut mask = Self::EMPTY;
        for element in iter {
            mask.insert(element);
        }
        mask
    }
}

/// A running minimum and maximum over values that may not all be recorded.
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub struct Extremes {
    min: Option<f32>,
    max: Option<f32>,
}

impl Extremes {
    /// Records a value, ignoring one that is not a number.
    pub fn observe(&mut self, value: f32) {
        if !value.is_finite() {
            return;
        }
        self.min = Some(match self.min {
            Some(current) => current.min(value),
            None => value,
        });
        self.max = Some(match self.max {
            Some(current) => current.max(value),
            None => value,
        });
    }

    /// The smallest value recorded.
    #[must_use]
    pub const fn min(self) -> Option<f32> {
        self.min
    }

    /// The largest value recorded.
    #[must_use]
    pub const fn max(self) -> Option<f32> {
        self.max
    }

    /// Returns true when no value in the range can satisfy `at_least`.
    ///
    /// This is what lets a filter on temperature factor skip a whole chunk.
    #[must_use]
    pub fn excludes_above(self, at_least: f32) -> bool {
        self.max.is_some_and(|max| max < at_least)
    }

    /// Returns true when no value in the range can satisfy `at_most`.
    #[must_use]
    pub fn excludes_below(self, at_most: f32) -> bool {
        self.min.is_some_and(|min| min > at_most)
    }
}

/// Everything a chunk knows about itself without reading a row.
#[derive(Clone, Copy, PartialEq, Debug, Default)]
pub struct AtomChunkStats {
    /// The box enclosing every recorded position.
    pub bounds: Aabb,
    /// The elements present.
    pub elements: ElementMask,
    /// The model every atom belongs to. Constant within a chunk by construction.
    pub model: u32,
    /// The lowest residue any atom belongs to.
    pub residue_min: u32,
    /// The highest residue any atom belongs to.
    pub residue_max: u32,
    /// The range of recorded temperature factors.
    pub b_factor: Extremes,
    /// The range of recorded occupancies.
    pub occupancy: Extremes,
    /// Whether any atom carries an alternate-location label.
    pub has_altloc: bool,
    /// Whether any atom's position was not recorded.
    pub has_missing_coords: bool,
    /// Whether any atom is hydrogen.
    pub has_hydrogen: bool,
}

impl AtomChunkStats {
    /// Returns true when no atom in the chunk can be `element`.
    #[must_use]
    pub const fn excludes_element(&self, element: Element) -> bool {
        !self.elements.contains(element)
    }

    /// Returns true when no atom in the chunk can be within `distance` of `region`.
    #[must_use]
    pub fn excludes_region(&self, region: &Aabb, distance: f32) -> bool {
        !self.bounds.within(region, distance)
    }

    /// Returns true when no atom in the chunk belongs to a residue in `range`.
    #[must_use]
    pub const fn excludes_residues(&self, first: u32, last: u32) -> bool {
        self.residue_max < first || self.residue_min > last
    }
}

#[cfg(test)]
#[path = "stats_tests.rs"]
mod tests;
