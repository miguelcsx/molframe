//! The windows that decide nucleic pucker class and suite sanity.
//!
//! A policy table is text under the shared table format. Each window header is
//! `name: low high` in degrees running upward, and each row is one staggered
//! gamma window `low high`:
//!
//! ```text
//! c3_endo_phase: -18 54
//! c2_endo_phase: 126 198
//! c3_endo_delta: 55 110
//! c2_endo_delta: 120 175
//! 30 90
//! ```
//!
//! Every window must be present and at least one gamma row is required; the
//! library compiles in no window.
//!
//! Complexity: decoding is O(rows) once, shared through an `Arc`; a window test
//! is O(1) and a gamma test is O(gamma windows).

use molframe_chem::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};

use super::nucleic_scoring::{PUCKER_C2_ENDO, PUCKER_C3_ENDO, PUCKER_OTHER};

/// An angular window in degrees from `low` to `high`, running upward.
pub type Window = [f64; 2];

const WHAT: &str = "nucleic policy table";

/// Windows that decide pucker class and suite sanity.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NucleicPolicy {
    /// Pseudorotation phase window of C3'-endo.
    pub c3_endo_phase: Window,
    /// Pseudorotation phase window of C2'-endo.
    pub c2_endo_phase: Window,
    /// Delta window consistent with a C3'-endo sugar.
    pub c3_endo_delta: Window,
    /// Delta window consistent with a C2'-endo sugar.
    pub c2_endo_delta: Window,
    /// Staggered gamma windows.
    pub gamma: Vec<Window>,
    seen: [bool; 4],
}

pub(super) fn within(angle: f64, window: Window) -> bool {
    (angle - window[0]).rem_euclid(360.0) <= (window[1] - window[0]).rem_euclid(360.0)
}

fn window(line: usize, value: &str) -> Result<Window, TableTextError> {
    let row = Row::new(WHAT, line);
    let [low, high] = row.words(value, "two numbers `low high`")?;
    Ok([row.finite(low)?, row.finite(high)?])
}

impl NucleicPolicy {
    /// A policy from its windows, for tables built in memory.
    #[must_use]
    pub const fn new(
        c3_endo_phase: Window,
        c2_endo_phase: Window,
        c3_endo_delta: Window,
        c2_endo_delta: Window,
        gamma: Vec<Window>,
    ) -> Self {
        Self {
            c3_endo_phase,
            c2_endo_phase,
            c3_endo_delta,
            c2_endo_delta,
            gamma,
            seen: [true; 4],
        }
    }

    pub(super) fn pucker_class(&self, phase: f64) -> u8 {
        if within(phase, self.c3_endo_phase) {
            PUCKER_C3_ENDO
        } else if within(phase, self.c2_endo_phase) {
            PUCKER_C2_ENDO
        } else {
            PUCKER_OTHER
        }
    }

    pub(super) fn delta_agrees(&self, class: u8, delta: f64) -> bool {
        match class {
            PUCKER_C3_ENDO => within(delta, self.c3_endo_delta),
            PUCKER_C2_ENDO => within(delta, self.c2_endo_delta),
            _ => within(delta, self.c3_endo_delta) || within(delta, self.c2_endo_delta),
        }
    }
}

impl TableLines for NucleicPolicy {
    fn header(&mut self, line: usize, key: &str, value: &str) -> Result<bool, TableTextError> {
        let (slot, index) = match key {
            "c3_endo_phase" => (&mut self.c3_endo_phase, 0),
            "c2_endo_phase" => (&mut self.c2_endo_phase, 1),
            "c3_endo_delta" => (&mut self.c3_endo_delta, 2),
            "c2_endo_delta" => (&mut self.c2_endo_delta, 3),
            _ => return Ok(false),
        };
        if self.seen[index] {
            return Err(TableTextError::at(
                WHAT,
                line,
                format!("`{key}` is given twice"),
            ));
        }
        *slot = window(line, value)?;
        self.seen[index] = true;
        Ok(true)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        self.gamma.push(window(line, text)?);
        Ok(())
    }
}

impl TableKind for NucleicPolicy {
    const KIND: &'static str = "nucleic";
    const EXTENSION: &'static str = "nucleic";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if !table.seen.iter().all(|seen| *seen) {
            return Err(TableTextError::at(
                WHAT,
                0,
                "c3_endo_phase, c2_endo_phase, c3_endo_delta and c2_endo_delta are all required",
            ));
        }
        if table.gamma.is_empty() {
            return Err(TableTextError::at(
                WHAT,
                0,
                "the table has no gamma window rows",
            ));
        }
        Ok(Provided { identity, table })
    }
}

#[cfg(test)]
#[path = "nucleic_policy_tests.rs"]
mod tests;
