//! Non-crystallographic operators of the `MTRIX1..3` records.

use crate::PdbHeaders;
use crate::fixed;

/// One non-crystallographic operator, from `MTRIX1..3`.
#[derive(Clone, Debug, PartialEq)]
pub struct Mtrix {
    /// Operator serial number.
    pub serial: u32,
    /// Row-major rotation.
    pub rotation: [[f64; 3]; 3],
    /// Translation in ångström.
    pub translation: [f64; 3],
    /// Whether the transformed copy is already present in the coordinates.
    pub given: bool,
}

/// Serial, the three rows read so far, and the given flag.
type PartialOperator = (u32, [Option<[f64; 4]>; 3], bool);

impl PdbHeaders {
    /// The `MTRIX` operators in deposited order.
    ///
    /// A serial that lacks any of its three rows is dropped rather than
    /// completed with zeros: a partial matrix would silently collapse atoms.
    #[must_use]
    pub fn mtrix(&self) -> Vec<Mtrix> {
        let mut rows: Vec<PartialOperator> = Vec::new();
        for (row, tag) in ["MTRIX1", "MTRIX2", "MTRIX3"].into_iter().enumerate() {
            for record in self.named(tag) {
                let line = record.line();
                let Some(serial) = fixed::integer(line, 8, 10).and_then(|n| u32::try_from(n).ok())
                else {
                    continue;
                };
                let (Some(a), Some(b), Some(c), Some(t)) = (
                    fixed::real(line, 11, 20),
                    fixed::real(line, 21, 30),
                    fixed::real(line, 31, 40),
                    fixed::real(line, 46, 55),
                ) else {
                    continue;
                };
                let given = fixed::text(line, 60, 60) == "1";
                let position = rows.iter().position(|entry| entry.0 == serial);
                let position = if let Some(position) = position {
                    position
                } else {
                    rows.push((serial, [None; 3], given));
                    rows.len() - 1
                };
                rows[position].1[row] = Some([a, b, c, t]);
            }
        }
        rows.into_iter()
            .filter_map(|(serial, rows, given)| {
                let [Some(x), Some(y), Some(z)] = rows else {
                    return None;
                };
                Some(Mtrix {
                    serial,
                    rotation: [[x[0], x[1], x[2]], [y[0], y[1], y[2]], [z[0], z[1], z[2]]],
                    translation: [x[3], y[3], z[3]],
                    given,
                })
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "mtrix_tests.rs"]
mod tests;
