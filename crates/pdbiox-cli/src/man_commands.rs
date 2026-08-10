//! Manual-page generation from the live declarative command tree.

use crate::exit::Exit;
use clap::Command;
use clap::CommandFactory as _;
use std::io::Write as _;
use std::path::Path;

pub(crate) fn completions(shell: crate::CompletionShell) -> Exit {
    let generator = match shell {
        crate::CompletionShell::Bash => clap_complete::Shell::Bash,
        crate::CompletionShell::Elvish => clap_complete::Shell::Elvish,
        crate::CompletionShell::Fish => clap_complete::Shell::Fish,
        crate::CompletionShell::PowerShell => clap_complete::Shell::PowerShell,
        crate::CompletionShell::Zsh => clap_complete::Shell::Zsh,
    };
    clap_complete::generate(
        generator,
        &mut crate::Cli::command(),
        "pdbiox",
        &mut std::io::stdout().lock(),
    );
    Exit::Success
}

pub(crate) fn generate(command: &Command, outdir: &Path) -> Exit {
    if !outdir.is_dir() {
        eprintln!("manual-page destination must be an existing directory");
        return Exit::Usage;
    }
    match render_tree(command, outdir, "pdbiox") {
        Ok(count) => {
            println!("generated {count} manual pages in {}", outdir.display());
            Exit::Success
        }
        Err(error) => {
            eprintln!("manual-page generation failed: {error}");
            Exit::Failure
        }
    }
}

fn render_tree(command: &Command, outdir: &Path, name: &str) -> std::io::Result<usize> {
    let path = outdir.join(format!("{name}.1"));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    clap_mangen::Man::new(command.clone()).render(&mut file)?;
    file.flush()?;
    let mut count = 1;
    for subcommand in command.get_subcommands() {
        let child_name = format!("{name}-{}", subcommand.get_name());
        count += render_tree(subcommand, outdir, &child_name)?;
    }
    Ok(count)
}
