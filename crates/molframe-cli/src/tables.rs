//! Reference tables for commands: one resolution, one refusal.
//!
//! Every command that needs a radius or reach table resolves it here, from the
//! command's own `--radii <NAME|PATH>` when given and from the profile
//! otherwise, against the process-wide store. A table that cannot be supplied is
//! refused with its `E7102` or `E2201` finding, which names how to supply it.
//!
//! Cost: one profile stat and, on first use, one table decode per command.

use crate::exit::Exit;
use molframe::chemistry::{
    CovalentRadii, Provided, RadiusSelection, RadiusSource, TableError, TableRequest, TableStore,
};
use molframe::{Diagnostic, Rendered};
use std::sync::Arc;

/// Prints the finding for `error` and returns the usage exit.
pub(crate) fn refuse(error: TableError) -> Exit {
    eprintln!("{}", Rendered::new(&Diagnostic::from(error)));
    Exit::Usage
}

/// The van der Waals radii a command measures with.
///
/// The one rule of [`RadiusSelection`]: `--radii`, else the table the configured policy's
/// `vdw_radii` names, else the profile; `--radii` naming another table than a configured
/// policy is the typed `E7103` refusal.
pub(crate) fn radius_source(
    argument: Option<&str>,
    context: &crate::report::Context,
) -> Result<RadiusSource, Exit> {
    RadiusSelection::choose(argument, context.radii_in_policy)
        .resolve(context.policy.vdw_radii, &TableStore::global())
        .map_err(refuse)
}

/// The covalent radii the profile selects.
pub(crate) fn covalent_radii() -> Result<Arc<Provided<CovalentRadii>>, Exit> {
    TableRequest::<CovalentRadii>::profile()
        .resolve(&TableStore::global(), "--covalent-radii")
        .map_err(refuse)
}
