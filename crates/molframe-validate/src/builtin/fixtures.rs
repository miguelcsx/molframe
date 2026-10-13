//! Synthetic reference tables for the tests, built from small closed-form rules.
//!
//! Nothing here is scientific data: each table only has the shape the decoders
//! require, with values chosen so a test can predict the outcome exactly.

use std::sync::Arc;

use molframe_chem::{Provided, TableIdentity, TableKind, TableRequest, TableStore, join_binary};

use super::ReferenceSet;
use super::nucleic_policy::NucleicPolicy;
use super::reference_data::{RAMA_BINS, RamaGrids, RotamerTables};
use super::reference_tables::{ChiDefinitions, GeometryTargets};

pub(super) fn identity(name: &str) -> TableIdentity {
    TableIdentity::supplied(name, "1", "synthetic test fixture", "CC0-1.0")
}

/// Six grids: every cell favoured except a block of forbidden cells around
/// phi = -179, psi = 1 in the general grid.
pub(super) fn rama_bytes() -> Vec<u8> {
    let cells = RAMA_BINS * RAMA_BINS;
    let mut payload = vec![6u8];
    for grid in 0..6 {
        payload.extend_from_slice(&[0u8; 8]);
        let mut body = vec![255u8; cells];
        if grid == 0 {
            let phi = 0;
            let psi = RAMA_BINS / 2;
            body[phi * RAMA_BINS + psi] = 0;
        }
        payload.extend_from_slice(&body);
    }
    join_binary(
        *b"MFR2",
        &identity("synthetic-rama"),
        "favoured: 0.02\noutlier: 0.0005 0.001 0.001 0.002 0.001 0.001\n",
        &payload,
    )
    .unwrap_or_default()
}

/// One rotamer table named `val` with three 120 degree bins: outlier, favoured, allowed.
pub(super) fn rotamer_bytes() -> Vec<u8> {
    let mut payload = vec![1u8, 3];
    payload.extend_from_slice(b"val");
    payload.push(1);
    payload.extend_from_slice(&[3, 0, 0]);
    payload.extend_from_slice(&3u32.to_le_bytes());
    payload.extend_from_slice(&[4, 6, 5]);
    join_binary(*b"MFT2", &identity("synthetic-rotamer"), "", &payload).unwrap_or_default()
}

pub(super) const GEOMETRY: &str = "name: synthetic-geometry\nversion: 1\n\
    source: synthetic test fixture\nlicence: CC0-1.0\n\
    bond\tALA\tN-CA\t1.46\t0.02\n\
    bond\tALA\tCA-C\t1.52\t0.02\n\
    angle\tALA\tN-CA-C\t111.0\t2.0\n\
    link\t*\tC-N\t1.33\t0.02\n\
    link\t*\tCA-C-N\t116.0\t2.0\n\
    link\tPRO\tC-N\t1.34\t0.02\n";

pub(super) const CHI: &str = "name: synthetic-chi\nversion: 1\n\
    source: synthetic test fixture\nlicence: CC0-1.0\n\
    VAL\tval\tN CA CB CG1\n";

pub(super) const NUCLEIC: &str = "name: synthetic-nucleic\nversion: 1\n\
    source: synthetic test fixture\nlicence: CC0-1.0\n\
    c3_endo_phase: -18 54\nc2_endo_phase: 126 198\n\
    c3_endo_delta: 55 110\nc2_endo_delta: 120 175\n\
    30 90\n150 210\n270 330\n";

fn given<T: TableKind>(bytes: &[u8]) -> TableRequest<T> {
    match T::decode(bytes, identity("fallback")) {
        Ok(provided) => TableRequest::given(provided),
        Err(_) => TableRequest::named("fixture-failed-to-decode"),
    }
}

/// A set with every table given in memory, resolved through an empty store.
pub(super) fn reference() -> ReferenceSet {
    let mut set = ReferenceSet::with_store(Arc::new(TableStore::new(
        std::env::temp_dir().join("molframe-no-such-tables"),
    )));
    set.rama = given::<RamaGrids>(&rama_bytes());
    set.rotamers = given::<RotamerTables>(&rotamer_bytes());
    set.geometry = given::<GeometryTargets>(GEOMETRY.as_bytes());
    set.chi = given::<ChiDefinitions>(CHI.as_bytes());
    set.nucleic = given::<NucleicPolicy>(NUCLEIC.as_bytes());
    set
}

/// The nucleic policy of [`NUCLEIC`] with `c3_endo_delta` replaced.
pub(super) fn nucleic_with_delta(delta: [f64; 2]) -> TableRequest<NucleicPolicy> {
    let policy = NucleicPolicy::new(
        [-18.0, 54.0],
        [126.0, 198.0],
        delta,
        [120.0, 175.0],
        vec![[30.0, 90.0], [150.0, 210.0], [270.0, 330.0]],
    );
    TableRequest::given(Provided {
        identity: identity("strict-nucleic"),
        table: policy,
    })
}
