//! Routes a parsed chemistry subcommand to the kernel that answers it.

use super::{peoe, smarts};
use crate::ChemCommand;
use crate::exit::Exit;
use crate::report::Context;

pub(crate) fn execute(command: ChemCommand, context: Context) -> Exit {
    match command {
        ChemCommand::Peoe { args } => peoe::peoe(&args, context),
        ChemCommand::Smarts { args } => smarts::smarts(&args, context),
    }
}
