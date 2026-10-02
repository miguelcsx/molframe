//! CI allocation ceilings for representative workflows.
//!
//! Every case runs in its own process because the counting allocator in the
//! binary is process-global, and one at a time: worker threads start lazily,
//! and a host busy with sibling processes starts them in a different order,
//! which moves the small counts by a few allocations. A ceiling is the largest
//! count observed over 200 sequential runs of the debug binary, times 1.05,
//! rounded up. Raising a ceiling means editing this table and justifying the
//! increase in the commit body.

use std::process::Command;
use std::sync::{Mutex, PoisonError};

/// Case name and the most allocations the case may perform.
const CEILINGS: &[(&str, u64)] = &[
    ("mmcif_read_medium", 16_997),
    ("bcif_read_medium", 24_602),
    ("selection_medium", 63),
    ("contacts_medium", 65),
    ("sasa_medium", 527),
    ("rmsd_medium", 38),
    ("dlpack_coordinates", 56),
    ("coordinate_handoff", 0),
];

/// Keeps the measured processes from overlapping across the tests below.
static ONE_AT_A_TIME: Mutex<()> = Mutex::new(());

struct Measurement {
    allocation_count: u64,
    result_digest: u64,
}

fn measure(case: &str) -> Measurement {
    let _turn = ONE_AT_A_TIME.lock().unwrap_or_else(PoisonError::into_inner);
    let output = Command::new(env!("CARGO_BIN_EXE_molframe-resource-bench"))
        .arg(case)
        .output()
        .expect("the resource bench binary starts");
    assert!(
        output.status.success(),
        "{case} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let record: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("the bench prints one JSON record");
    let field = |name: &str| {
        record[name]
            .as_u64()
            .unwrap_or_else(|| panic!("{case}: field {name} is not an unsigned integer"))
    };
    Measurement {
        allocation_count: field("allocation_count"),
        result_digest: field("result_digest"),
    }
}

/// The fewest allocations a case performs over three runs.
///
/// The counter is process-wide, so a pool thread finishing its own start-up
/// inside the measured window adds to the count; it never subtracts. The
/// case itself is deterministic, so the minimum of a few runs is its count.
fn allocations(case: &str) -> u64 {
    (0..3)
        .map(|_| measure(case).allocation_count)
        .min()
        .unwrap_or(u64::MAX)
}

#[test]
fn allocation_counts_stay_under_their_ceilings() {
    for &(case, ceiling) in CEILINGS {
        let measured = allocations(case);
        assert!(
            measured <= ceiling,
            "{case} allocated {measured} times, above its ceiling of {ceiling}"
        );
    }
}

#[test]
fn coordinate_handoff_allocates_nothing() {
    assert_eq!(allocations("coordinate_handoff"), 0);
}

#[test]
fn results_are_digest_stable() {
    for &(case, _) in CEILINGS {
        let first = measure(case).result_digest;
        let second = measure(case).result_digest;
        assert_eq!(first, second, "{case} produced different result digests");
    }
}
