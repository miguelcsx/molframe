//! The metal ions of structural biology: which elements they are and how big.
//!
//! Surface and steric clash work meet the same ions and need two radii for
//! them: a van der Waals radius that sizes the ion as a solvent obstacle, and
//! an ionic radius that sizes a coordinated cation as a steric contact. Which
//! elements are metals, the oxidation state and coordination each most often has
//! in deposited macromolecules, and their van der Waals radius are a table like
//! any other, selected by name (profile line `metals`, files `<name>.metals`),
//! never compiled in:
//!
//! ```text
//! # element typical-charge typical-coordination spin(-|HS|LS) vdw-radius
//! Fe 2 VI HS 2.44
//! ```
//!
//! The ionic radius of a site is read from the ionic radius table for exactly
//! that charge, coordination and spin, so the two tables agree by lookup, not by
//! copying. A metal with no such entry has no ionic radius.
//!
//! Cost: decoding is O(rows); a lookup is one array index.

use crate::{
    IonicRadii, IonicSpin, Provided, Row, TableError, TableIdentity, TableKind, TableLines,
    TableRequest, TableStore, TableTextError, read_table, text_of,
};
use molframe_core::Element;
use std::sync::Arc;

const WHAT: &str = "metal site table";

/// Which radius of a metal ion a consumer wants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MetalRadiusKind {
    /// The van der Waals radius, the size of the ion as a solvent obstacle.
    Vdw,
    /// The ionic radius at the metal's typical charge and coordination, the
    /// size of a coordinated cation as a steric contact.
    Ionic,
}

/// The oxidation state and coordination a metal usually has in macromolecules.
#[derive(Clone, Debug, PartialEq)]
pub struct MetalSite {
    /// Typical formal charge.
    pub charge: i8,
    /// Typical coordination in Shannon's notation.
    pub coordination: Box<str>,
    /// Spin state the typical site takes.
    pub spin: IonicSpin,
    /// Van der Waals radius in ångström.
    pub vdw_radius: f32,
}

/// A validated table of metal sites, one per element.
#[derive(Clone, Debug, PartialEq)]
pub struct MetalSites {
    sites: Vec<Option<MetalSite>>,
}

impl Default for MetalSites {
    fn default() -> Self {
        Self {
            sites: vec![None; 119],
        }
    }
}

impl TableKind for MetalSites {
    const KIND: &'static str = "metals";
    const EXTENSION: &'static str = "metals";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if table.sites.iter().all(Option::is_none) {
            return Err(TableTextError::at(WHAT, 0, "the table lists no metal"));
        }
        Ok(Provided { identity, table })
    }
}

impl TableLines for MetalSites {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let row = Row::new(WHAT, line);
        let [symbol, charge, coordination, spin, vdw] =
            row.words(text, "`element charge coordination spin vdw-radius`")?;
        let element = row.element(symbol)?;
        let site = MetalSite {
            charge: charge
                .parse()
                .map_err(|_| row.fail(format!("`{charge}` is not a charge")))?,
            coordination: coordination.into(),
            spin: IonicSpin::from_word(spin).map_err(|message| row.fail(message))?,
            vdw_radius: row.bounded::<f32>(vdw, 20.0, "the van der Waals radius")?,
        };
        row.claim(
            self.sites.get_mut(element.atomic_number() as usize),
            symbol,
            site,
        )
    }
}

impl MetalSites {
    /// The typical site of a metal, or `None` for any other element.
    #[must_use]
    pub fn site(&self, element: Element) -> Option<&MetalSite> {
        self.sites
            .get(element.atomic_number() as usize)
            .and_then(Option::as_ref)
    }

    /// Whether the table lists `element` as a metal.
    #[must_use]
    pub fn is_metal(&self, element: Element) -> bool {
        self.site(element).is_some()
    }

    /// The elements listed, in atomic-number order.
    pub fn elements(&self) -> impl Iterator<Item = Element> + '_ {
        self.sites
            .iter()
            .enumerate()
            .filter(|(_, site)| site.is_some())
            .filter_map(|(z, _)| u8::try_from(z).ok().map(Element::from_atomic_number))
    }

    /// A radius of a metal ion in ångström, or `None` for an element the table
    /// does not list or that has no value of that kind. The ionic radius needs
    /// `ionic`; without it none exists.
    #[must_use]
    pub fn radius(
        &self,
        element: Element,
        kind: MetalRadiusKind,
        ionic: Option<&IonicRadii>,
    ) -> Option<f32> {
        let site = self.site(element)?;
        match kind {
            MetalRadiusKind::Vdw => Some(site.vdw_radius),
            MetalRadiusKind::Ionic => ionic?
                .entries(element)
                .iter()
                .find(|entry| {
                    entry.charge == site.charge
                        && entry.coordination == site.coordination
                        && entry.spin == site.spin
                })
                .map(|entry| entry.ionic_radius),
        }
    }
}

/// The resolved tables a metal fallback draws on.
#[derive(Clone, Debug, PartialEq)]
pub struct MetalTables {
    /// Which elements are metals and their van der Waals radii.
    pub sites: Arc<Provided<MetalSites>>,
    /// The ionic radius catalogue, present when ionic radii are asked for.
    pub ionic: Option<Arc<Provided<IonicRadii>>>,
}

impl MetalTables {
    /// Resolves the metal table, and the ionic catalogue when `with_ionic`.
    ///
    /// # Errors
    ///
    /// [`TableError`] when a table is not selected or is refused.
    pub fn resolve(
        sites: &TableRequest<MetalSites>,
        ionic: Option<&TableRequest<IonicRadii>>,
        store: &TableStore,
    ) -> Result<Self, TableError> {
        Ok(Self {
            sites: sites.resolve(store, "--metals")?,
            ionic: ionic
                .map(|request| request.resolve(store, "--ionic-radii"))
                .transpose()?,
        })
    }

    /// A radius of `element`; see [`MetalSites::radius`].
    #[must_use]
    pub fn radius(&self, element: Element, kind: MetalRadiusKind) -> Option<f32> {
        let ionic = self.ionic.as_ref().map(|provided| &provided.table);
        self.sites.table.radius(element, kind, ionic)
    }

    /// Identities of the tables in force.
    #[must_use]
    pub fn identities(&self) -> Vec<TableIdentity> {
        let mut all = vec![self.sites.identity.clone()];
        all.extend(self.ionic.iter().map(|ionic| ionic.identity.clone()));
        all
    }
}

#[cfg(test)]
#[path = "metal_tests.rs"]
mod tests;
