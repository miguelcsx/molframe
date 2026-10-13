//! Struct-of-arrays form of a pair list, ready to hand to column consumers.

use super::NeighborPair;
use crate::CrossPair;

/// A pair list split into index columns and Euclidean distances.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PairColumns {
    /// Lower atom index of every pair.
    pub first: Vec<u32>,
    /// Higher atom index of every pair.
    pub second: Vec<u32>,
    /// Euclidean separation of every pair, in ångström.
    pub distance: Vec<f32>,
}

impl PairColumns {
    /// Splits `pairs` into columns, taking the square root of each distance.
    #[must_use]
    pub fn from_pairs(pairs: &[NeighborPair]) -> Self {
        let mut columns = Self {
            first: Vec::with_capacity(pairs.len()),
            second: Vec::with_capacity(pairs.len()),
            distance: Vec::with_capacity(pairs.len()),
        };
        for pair in pairs {
            columns.first.push(pair.first);
            columns.second.push(pair.second);
            columns.distance.push(pair.distance_squared.sqrt());
        }
        columns
    }

    /// Splits cross-set `pairs` into columns; indices stay local to each set.
    #[must_use]
    pub fn from_cross(pairs: &[CrossPair]) -> Self {
        let mut columns = Self {
            first: Vec::with_capacity(pairs.len()),
            second: Vec::with_capacity(pairs.len()),
            distance: Vec::with_capacity(pairs.len()),
        };
        for pair in pairs {
            columns.first.push(pair.first);
            columns.second.push(pair.second);
            columns.distance.push(pair.distance_squared.sqrt());
        }
        columns
    }
}
