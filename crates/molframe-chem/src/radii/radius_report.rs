//! The radius policy a calculation actually applied, stated in full.
//!
//! A result that says only "Bondi radii" cannot be reproduced against another
//! release or compared with a run that used the caller's own table. This
//! report names every table that fed the radii, with its version, source and
//! licence, together with the choices that decided which atoms took part.

use crate::{HeteroPolicy, MetalFallback, MetalTables, RadiusPolicy, TableIdentity};
use core::fmt;

/// What a [`RadiusPolicy`] resolves to, table by table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AppliedRadiusPolicy {
    /// The table that supplied the radii.
    pub radii: TableIdentity,
    /// Which fallback covered metals the table lacks.
    pub metals: MetalFallback,
    /// The tables behind that fallback, empty when it is `none`.
    pub metal_tables: Vec<TableIdentity>,
    /// Which heteroatom classes took part.
    pub hetero: HeteroPolicy,
}

impl RadiusPolicy {
    /// The policy with every table it resolves to named.
    #[must_use]
    pub fn applied(&self) -> AppliedRadiusPolicy {
        AppliedRadiusPolicy {
            radii: self.source.identity(),
            metals: self.metals,
            metal_tables: self
                .metal_tables
                .as_ref()
                .map_or_else(Vec::new, MetalTables::identities),
            hetero: self.hetero,
        }
    }
}

impl fmt::Display for AppliedRadiusPolicy {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(out, "radii {}; metals {}", self.radii, self.metals)?;
        if !self.metal_tables.is_empty() {
            let tables: Vec<String> = self.metal_tables.iter().map(ToString::to_string).collect();
            write!(out, " [{}]", tables.join(", "))?;
        }
        write!(out, "; hetero {}", self.hetero)
    }
}

#[cfg(test)]
#[path = "radius_report_tests.rs"]
mod tests;
