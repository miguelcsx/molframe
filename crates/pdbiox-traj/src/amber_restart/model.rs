//! Source metadata retained with a canonical trajectory frame.

use crate::{AmberRestartLayout, Timestep};

/// One formatted AMBER restart/inpcrd record.
#[derive(Clone, Debug, PartialEq)]
pub struct AmberRestart {
    /// Title line, preserved without its line terminator.
    pub title: Box<str>,
    /// Explicit on-disk field layout.
    pub layout: AmberRestartLayout,
    /// Values in canonical pdbiox units: ångström, picosecond, and ångström/ps.
    pub timestep: Timestep,
}
