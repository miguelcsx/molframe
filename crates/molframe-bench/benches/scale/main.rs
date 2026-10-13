//! Declarative scale sweep: every case, over size tiers, at chosen worker counts.
//!
//! ```text
//! cargo bench -p molframe-bench --bench scale -- \
//!     [--filter text] [--max-tier large] [--workers 1,4] \
//!     [--budget-secs 30] [--out report.tsv]
//! ```
//!
//! Each row records the best of several repetitions, peak live heap above the
//! baseline, allocation count and a bit-exact result checksum. Rows for one
//! case and size at different worker counts must carry the same checksum; the
//! report flags any that do not. After the rows, one line per case and worker
//! count gives the fitted scaling exponent.

mod alloc;
mod compare_cases;
mod geom_cases;
mod harness;
mod interop_cases;
mod registry;
mod report;
mod spatial_cases;
mod traj_cases;

use harness::{Case, Settings};
use molframe_bench::Tier;

#[global_allocator]
static ALLOCATOR: alloc::CountingAllocator = alloc::CountingAllocator;

fn main() {
    let settings = match Settings::from_args(std::env::args().skip(1)) {
        Ok(settings) => settings,
        Err(message) => {
            eprintln!("scale: {message}");
            std::process::exit(2);
        }
    };
    let all: Vec<Case> = registry::all();
    let selected: Vec<&Case> = all
        .iter()
        .filter(|case| {
            settings
                .filter
                .as_deref()
                .is_none_or(|f| case.id.contains(f))
        })
        .collect();
    let rows = harness::run(&selected, &settings);
    report::emit(&rows, &settings, Tier::ALL);
}
