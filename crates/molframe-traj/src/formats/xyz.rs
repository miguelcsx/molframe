//! Reading the XYZ coordinate format.
//!
//! An XYZ frame is an atom count, a free-text comment, then one line per atom
//! giving its element symbol and Cartesian position. A file may hold several
//! frames back to back, which is how XYZ carries a trajectory. This reads them
//! all, and converting a frame to a bare coordinate [`Frame`] drops the elements
//! for the coordinate store.
//!
//! Runs in `O(atoms)` per frame.

use molframe_core::element::Element;
use std::fmt::Write;

use crate::trajectory::Frame;

/// One atom of an XYZ frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct XyzAtom {
    /// The element, or [`Element::UNKNOWN`] for an unrecognised symbol.
    pub element: Element,
    /// The atom position.
    pub position: [f32; 3],
}

/// One XYZ frame: its comment line and its atoms.
#[derive(Clone, Debug, PartialEq)]
pub struct XyzFrame {
    /// The free-text comment line.
    pub comment: String,
    /// The atoms in file order.
    pub atoms: Vec<XyzAtom>,
}

impl XyzFrame {
    /// The frame's positions as a coordinate [`Frame`], dropping the elements.
    #[must_use]
    pub fn to_frame(&self) -> Frame {
        Frame {
            positions: self.atoms.iter().map(|atom| atom.position).collect(),
        }
    }
}

/// Parses every frame of an XYZ file, returning `None` on a malformed record.
#[must_use]
pub fn parse_xyz(text: &str) -> Option<Vec<XyzFrame>> {
    let mut frames = Vec::new();
    let mut lines = text.lines();
    loop {
        let count_line = loop {
            match lines.next() {
                Some(line) if line.trim().is_empty() => {}
                Some(line) => break line,
                None => return Some(frames),
            }
        };
        let atom_count: usize = count_line.trim().parse().ok()?;
        let comment = lines.next()?.to_string();
        let mut atoms = Vec::with_capacity(atom_count);
        for _ in 0..atom_count {
            let mut fields = lines.next()?.split_whitespace();
            let symbol = fields.next()?;
            let x: f32 = fields.next()?.parse().ok()?;
            let y: f32 = fields.next()?.parse().ok()?;
            let z: f32 = fields.next()?.parse().ok()?;
            let element = match Element::from_symbol(symbol) {
                Some(element) => element,
                None => Element::UNKNOWN,
            };
            atoms.push(XyzAtom {
                element,
                position: [x, y, z],
            });
        }
        frames.push(XyzFrame { comment, atoms });
    }
}

/// Writes XYZ frames without dropping their elements or comments.
#[must_use]
pub fn write_xyz(frames: &[XyzFrame]) -> String {
    let mut output = String::new();
    for frame in frames {
        let _ = writeln!(output, "{}", frame.atoms.len());
        let _ = writeln!(output, "{}", frame.comment);
        for atom in &frame.atoms {
            let _ = writeln!(
                output,
                "{} {:.6} {:.6} {:.6}",
                atom.element.symbol(),
                atom.position[0],
                atom.position[1],
                atom.position[2],
            );
        }
    }
    output
}

#[cfg(test)]
#[path = "xyz_tests.rs"]
mod tests;
