//! Which van der Waals table a calculation uses when a policy is in play.
//!
//! The policy's `vdw_radii` names a radius set, and a calculation may also be told a
//! table directly (`radii`). Those are two answers to one question, so the rule that
//! combines them is stated once, here, and every consumer applies it:
//!
//! - a table asked for directly is the table used;
//! - with none, the table the policy names, found by that name through the
//!   [`TableStore`] and never by a name this crate knows;
//! - both given and naming different tables is [`TableError::Conflict`], never a
//!   silent choice of one;
//! - a policy-named table the store cannot supply is [`TableError::Missing`], never
//!   a quiet fall back to the profile's table.
//!
//! A caller with no policy of its own leaves the policy out, and the profile decides
//! as it always did.
//!
//! Cost: one store lookup, as resolving any request.

use crate::{
    RadiusSource, TableError, TableIdentity, TableKind, TableRequest, TableStore, VdwRadii,
};
use molframe_core::contract::RadiiSet;

/// How a calculation picks its radius table relative to a policy.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiusSelection {
    /// The table the policy names.
    Policy,
    /// This table; no policy was supplied, so none is consulted.
    Table(TableRequest<VdwRadii>),
    /// This table, which must be the one the policy names.
    Agreeing(TableRequest<VdwRadii>),
}

impl RadiusSelection {
    /// The selection a front end's optional `radii` argument makes, given whether the
    /// caller supplied a policy: the argument alone without one, the argument
    /// checked against the policy with one, and the policy's table when there is no
    /// argument. Neither argument nor policy leaves the profile to decide.
    #[must_use]
    pub fn choose(argument: Option<&str>, policy_supplied: bool) -> Self {
        match (argument, policy_supplied) {
            (None, true) => Self::Policy,
            (None, false) => Self::Table(TableRequest::Profile),
            (Some(name), false) => Self::Table(TableRequest::named(name)),
            (Some(name), true) => Self::Agreeing(TableRequest::named(name)),
        }
    }

    /// The table this selection applies under a policy that names `set`.
    ///
    /// # Errors
    ///
    /// [`TableError::Missing`] when nothing supplies the table the policy names,
    /// [`TableError::Conflict`] when the argument and the policy name different
    /// tables, or the error of reading or validating a file.
    pub fn resolve(&self, set: RadiiSet, store: &TableStore) -> Result<RadiusSource, TableError> {
        let source = match self {
            Self::Policy => RadiusSource::resolve(&TableRequest::named(set.name()), store)?,
            Self::Table(request) | Self::Agreeing(request) => {
                RadiusSource::resolve(request, store)?
            }
        };
        self.verify(&source.identity(), set)?;
        Ok(source)
    }

    /// Checks that `applied`, the table a calculation was handed, is the one this
    /// selection allows under a policy naming `set`.
    ///
    /// # Errors
    ///
    /// [`TableError::Conflict`] when the selection ties the table to the policy and
    /// the table is another.
    pub fn verify(&self, applied: &TableIdentity, set: RadiiSet) -> Result<(), TableError> {
        match self {
            Self::Table(_) => Ok(()),
            Self::Policy | Self::Agreeing(_) if applied.name == set.name() => Ok(()),
            Self::Policy | Self::Agreeing(_) => Err(TableError::Conflict {
                kind: VdwRadii::KIND,
                requested: applied.name.to_string(),
                policy: set.name().to_owned(),
            }),
        }
    }

    /// The table this selection applies under `set`, as a request that selects
    /// exactly it, for option sets that take a request.
    ///
    /// # Errors
    ///
    /// As [`RadiusSelection::resolve`].
    pub fn request(
        &self,
        set: RadiiSet,
        store: &TableStore,
    ) -> Result<TableRequest<VdwRadii>, TableError> {
        self.resolve(set, store).map(|source| source.request())
    }
}

#[cfg(test)]
#[path = "radius_selection_tests.rs"]
mod tests;
