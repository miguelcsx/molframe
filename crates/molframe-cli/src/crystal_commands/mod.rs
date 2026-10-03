//! Crystallographic data commands: reflections, density maps and symmetry mates.

mod maps;
mod mates;
mod mtz;

use crate::CrystalCommand;
use crate::exit::Exit;
use crate::report::Context;

pub(crate) fn execute(command: CrystalCommand, context: Context) -> Exit {
    match command {
        CrystalCommand::MtzInfo { input, table } => mtz::info(&input, table, context),
        CrystalCommand::MapStats { input, histogram } => {
            maps::statistics(&input, &histogram, context)
        }
        CrystalCommand::MapCorrelation {
            observed,
            calculated,
        } => maps::correlation(&observed, &calculated, context),
        CrystalCommand::Mates { args } => mates::mates(&args, context),
    }
}
