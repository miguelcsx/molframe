//! CI allocation ceilings for representative workflows.
//!
//! Every case runs in its own process because the counting allocator in the
//! binary is process-global. A ceiling is the largest count observed over 100
//! runs of the debug test binary, times 1.05, rounded up. The small cases
//! vary by tens of allocations from worker-thread start-up, so the maximum
//! rather than the median is the baseline. Raising a ceiling means editing
//! this table and justifying the increase in the commit body.

use std::process::Command;

/// Case name and the most allocations the case may perform.
const CEILINGS: &[(&str, u64)] = &[
    ("mmcif_read_medium", 16_997),
    ("bcif_read_medium", 24_701),
    ("selection_medium", 63),
    ("contacts_medium", 65),
    ("sasa_medium", 220),
    ("rmsd_medium", 38),
    ("dlpack_coordinates", 56),
    ("coordinate_handoff", 0),
];

struct Measurement {
    allocation_count: u64,
    result_digest: u64,
}

fn measure(case: &str) -> Measurement {
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

#[test]
fn allocation_counts_stay_under_their_ceilings() {
    for &(case, ceiling) in CEILINGS {
        let measured = measure(case).allocation_count;
        assert!(
            measured <= ceiling,
            "{case} allocated {measured} times, above its ceiling of {ceiling}"
        );
    }
}

#[test]
fn coordinate_handoff_allocates_nothing() {
    assert_eq!(measure("coordinate_handoff").allocation_count, 0);
}

#[test]
fn results_are_digest_stable() {
    for &(case, _) in CEILINGS {
        let first = measure(case).result_digest;
        let second = measure(case).result_digest;
        assert_eq!(first, second, "{case} produced different result digests");
    }
}
