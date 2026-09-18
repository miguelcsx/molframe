//! CLI projection over verified external-resource retrieval.

use crate::exit::Exit;
use crate::report::{Context, Json};
use molframe_adapters::VerifiedDownload;
use std::io::Write as _;
use std::path::Path;

const ID_PLACEHOLDER: &str = "{id}";

#[derive(Clone, Copy)]
pub(crate) struct NetworkOptions {
    pub(crate) max_bytes: u64,
    pub(crate) timeout_seconds: u64,
    pub(crate) redirect_limit: usize,
}

pub(crate) fn fetch(
    id: &str,
    url_template: &str,
    sha256: &str,
    options: NetworkOptions,
    output: &Path,
    context: Context,
) -> Exit {
    if !valid_id(id) || !url_template.contains(ID_PLACEHOLDER) {
        eprintln!("fetch requires an alphanumeric ID and --url-template containing {{id}}");
        return Exit::Usage;
    }
    let url = url_template.replace(ID_PLACEHOLDER, id);
    let download = match verified_download(&url, sha256, options) {
        Ok(download) => download,
        Err(error) => {
            eprintln!("fetch failed: {error}");
            return Exit::Input;
        }
    };
    let name = output.file_name().and_then(std::ffi::OsStr::to_str);
    let options = molframe::ReadOptions::new()
        .mode(context.mode)
        .missing_element_policy(context.missing_element_policy)
        .ambiguous_residue_boundary_policy(context.residue_boundary_policy);
    match molframe::read_bytes(download.bytes.clone(), name, &options) {
        Ok((_, findings)) => context.findings(&findings, &url),
        Err(findings) => {
            context.findings(&findings, &url);
            return Exit::of(&findings);
        }
    }
    write_download(output, &download, false, context)
}

pub(crate) fn update_ccd(
    url: &str,
    sha256: &str,
    options: NetworkOptions,
    version: &str,
    output: &Path,
    replace: bool,
    context: Context,
) -> Exit {
    if version.is_empty() {
        eprintln!("CCD update requires a non-empty explicit version");
        return Exit::Usage;
    }
    let download = match verified_download(url, sha256, options) {
        Ok(download) => download,
        Err(error) => {
            eprintln!("CCD download failed: {error}");
            return Exit::Input;
        }
    };
    // The downloaded bytes are validated before they are written, so this reads
    // from the buffer rather than from a path the facade could open.
    let input = molframe::InputBuffer::from_bytes(download.bytes.clone());
    match molframe::chem::read_ccd(&input, molframe::DictionaryVersion::new(version)) {
        Ok((_, findings)) => context.findings(&findings, url),
        Err(findings) => {
            context.findings(&findings, url);
            return Exit::of(&findings);
        }
    }
    write_download(output, &download, replace, context)
}

fn verified_download(
    url: &str,
    sha256: &str,
    options: NetworkOptions,
) -> Result<VerifiedDownload, molframe_adapters::DownloadError> {
    molframe_adapters::fetch_verified(
        url,
        sha256,
        molframe_adapters::DownloadOptions {
            max_bytes: options.max_bytes,
            timeout: std::time::Duration::from_secs(options.timeout_seconds),
            redirect_limit: options.redirect_limit,
        },
    )
}

fn write_download(
    output: &Path,
    download: &VerifiedDownload,
    replace: bool,
    context: Context,
) -> Exit {
    let result = if replace {
        atomic_replace(output, &download.bytes)
    } else {
        write_new(output, &download.bytes)
    };
    match result {
        Ok(()) => {
            if context.is_json() {
                let mut object = Json::new();
                object
                    .text("output", &output.display().to_string())
                    .number("bytes", download.bytes.len())
                    .text("sha256", &download.sha256);
                context.result(&object.finish());
            } else {
                context.result(&format!(
                    "wrote {} verified bytes to {} ({})",
                    download.bytes.len(),
                    output.display(),
                    download.sha256
                ));
            }
            Exit::Success
        }
        Err(error) => {
            eprintln!("could not write {}: {error}", output.display());
            Exit::Failure
        }
    }
}

fn write_new(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    file.write_all(bytes)?;
    file.sync_all()
}

fn atomic_replace(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty());
    let directory = match parent {
        Some(parent) => parent,
        None => Path::new("."),
    };
    let mut temporary = tempfile::NamedTempFile::new_in(directory)?;
    temporary.write_all(bytes)?;
    temporary.as_file().sync_all()?;
    temporary.persist(path).map_err(|error| error.error)?;
    Ok(())
}

fn valid_id(id: &str) -> bool {
    !id.is_empty()
        && id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_'))
}
