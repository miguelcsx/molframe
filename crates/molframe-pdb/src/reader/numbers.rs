//! Per-row numeric fields, decoded ahead of the sequential read.
//!
//! The text-to-number work of a coordinate row depends on nothing but the row,
//! so a window of rows decodes on the shared pool while everything that carries
//! state from row to row (chains, residues, interning, findings) stays a single
//! ordered walk that consumes the results.

use super::lines::Line;
use super::state::{charge_in, serial_in};
use crate::fixed;
use molframe_core::column::Presence;
use molframe_core::element::Element;
use molframe_core::execution::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

/// Rows decoded per scheduling block.
const BLOCK_LINES: usize = 2048;

/// The row-local numbers of one coordinate record.
#[derive(Clone, Copy)]
pub(super) struct AtomNumbers {
    pub(super) position: [Option<f64>; 3],
    pub(super) occupancy: Option<f64>,
    pub(super) b_factor: Option<f64>,
    pub(super) serial: u32,
    pub(super) charge: (i8, Presence),
    pub(super) declared: Option<Element>,
}

impl AtomNumbers {
    fn parse(text: &str) -> Self {
        Self {
            position: [
                fixed::real(text, 31, 38),
                fixed::real(text, 39, 46),
                fixed::real(text, 47, 54),
            ],
            occupancy: fixed::real(text, 55, 60),
            b_factor: fixed::real(text, 61, 66),
            serial: serial_in(text),
            charge: charge_in(text),
            declared: Element::from_symbol(fixed::text(text, 77, 78)),
        }
    }
}

/// Decodes the numbers of every coordinate row in `lines`, in row order.
///
/// A row that is not a coordinate record yields `None`, and so does every row
/// if a worker fails: the sequential read then decodes the row itself, which
/// is what it did before this window existed.
pub(super) fn decode_window(
    lines: &[Line<'_>],
    context: &ExecutionContext,
) -> Vec<Option<AtomNumbers>> {
    let plan = BlockPlan::new(lines.len(), BLOCK_LINES);
    let blocks = map_blocks_in(plan, context, |_, range| {
        lines
            .get(range)
            .unwrap_or_default()
            .iter()
            .map(|line| {
                matches!(fixed::record(line.text), "ATOM" | "HETATM")
                    .then(|| AtomNumbers::parse(line.text))
            })
            .collect::<Vec<_>>()
    });
    match blocks {
        Ok(blocks) => blocks.into_iter().flatten().collect(),
        Err(_) => Vec::new(),
    }
}
