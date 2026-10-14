//! CLI projection over the manifest-driven table cache.
//!
//! The commands parse arguments and print; the manifest, conversion, download
//! and validation are `molframe::data`. `list` and `verify` read local files
//! only; `fetch` is the one command that touches the network.

use crate::DataCommand;
use crate::exit::Exit;
use crate::report::{Context, Json, json_array};
use molframe::chemistry::TableStore;
use molframe::data::{DataError, Manifest, Status, Verified};
use std::path::{Path, PathBuf};

pub(crate) fn execute(command: DataCommand, context: Context) -> Exit {
    match command {
        DataCommand::List { manifest, dir } => with_manifest(&manifest, dir, |manifest, dir| {
            show(&molframe::data::list(manifest, dir), context);
            Exit::Success
        }),
        DataCommand::Verify {
            name,
            manifest,
            dir,
        } => with_manifest(
            &manifest,
            dir,
            |manifest, dir| match molframe::data::verify(manifest, dir, name.as_deref()) {
                Ok(rows) => {
                    show(&rows, context);
                    if rows.iter().all(|row| row.status == Status::Current) {
                        Exit::Success
                    } else {
                        Exit::Input
                    }
                }
                Err(error) => failed(&error),
            },
        ),
        DataCommand::Fetch {
            name,
            manifest,
            dir,
            replace,
            max_bytes,
            timeout_seconds,
            redirect_limit,
        } => with_manifest(&manifest, dir, |manifest, dir| {
            let options = molframe::adapters::DownloadOptions {
                max_bytes,
                timeout: std::time::Duration::from_secs(timeout_seconds),
                redirect_limit,
            };
            match molframe::data::fetch(manifest, &name, dir, replace, options) {
                Ok(stored) => {
                    context.result(&format!(
                        "stored {} ({} rows, {} bytes)",
                        stored.path.display(),
                        stored.rows,
                        stored.bytes
                    ));
                    Exit::Success
                }
                Err(error) => failed(&error),
            }
        }),
    }
}

fn with_manifest(
    manifest: &Path,
    dir: Option<PathBuf>,
    run: impl FnOnce(&Manifest, &Path) -> Exit,
) -> Exit {
    let manifest = match Manifest::load(manifest) {
        Ok(manifest) => manifest,
        Err(error) => return failed(&error),
    };
    match dir {
        Some(dir) => run(&manifest, &dir),
        None => run(&manifest, TableStore::global().directory()),
    }
}

fn failed(error: &DataError) -> Exit {
    eprintln!("{error}");
    match error {
        DataError::NoEntry { .. } => Exit::Usage,
        DataError::Manifest { .. } | DataError::Convert { .. } => Exit::Schema,
        DataError::Read { .. } | DataError::Download { .. } => Exit::Input,
        DataError::Exists { .. } | DataError::Write { .. } => Exit::Failure,
    }
}

fn show(rows: &[Verified], context: Context) {
    if context.is_json() {
        let objects: Vec<String> = rows
            .iter()
            .map(|row| {
                let mut object = Json::new();
                object
                    .text("name", &row.entry.name)
                    .text("version", &row.entry.version)
                    .text("licence", &row.entry.licence)
                    .text("path", &row.path.display().to_string())
                    .text("status", &row.status.to_string());
                object.finish()
            })
            .collect();
        context.result(&json_array(&objects));
    } else {
        for row in rows {
            context.result(&format!(
                "{}\t{}\t{}\t{}",
                row.entry.name, row.entry.version, row.entry.licence, row.status
            ));
        }
    }
}
