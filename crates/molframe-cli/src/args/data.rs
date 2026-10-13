//! Arguments for the manifest-driven reference-table cache.

use clap::Subcommand;
use std::path::PathBuf;

#[derive(Subcommand, Debug)]
pub(crate) enum DataCommand {
    /// List the entries of a manifest and whether each is in the cache.
    List {
        /// TOML manifest naming the entries.
        #[arg(long)]
        manifest: PathBuf,
        /// Cache directory; defaults to `MOLFRAME_DATA_DIR` or the user cache.
        #[arg(long)]
        dir: Option<PathBuf>,
    },
    /// Download one entry, verify its pinned SHA-256, convert and store it.
    Fetch {
        /// Entry name in the manifest.
        name: String,
        /// TOML manifest naming the entries.
        #[arg(long)]
        manifest: PathBuf,
        /// Cache directory; defaults to `MOLFRAME_DATA_DIR` or the user cache.
        #[arg(long)]
        dir: Option<PathBuf>,
        /// Replace an existing cache entry.
        #[arg(long)]
        replace: bool,
        /// Maximum accepted download size in bytes.
        #[arg(long, default_value_t = 268_435_456)]
        max_bytes: u64,
        /// Whole-request timeout in seconds.
        #[arg(long, default_value_t = 120)]
        timeout_seconds: u64,
        /// Maximum HTTP redirects; zero refuses redirects.
        #[arg(long, default_value_t = 5)]
        redirect_limit: usize,
    },
    /// Check cached entries against the manifest without any network access.
    Verify {
        /// Entry name; every entry when omitted.
        name: Option<String>,
        /// TOML manifest naming the entries.
        #[arg(long)]
        manifest: PathBuf,
        /// Cache directory; defaults to `MOLFRAME_DATA_DIR` or the user cache.
        #[arg(long)]
        dir: Option<PathBuf>,
    },
}
