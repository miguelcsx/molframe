//! Per-residue nucleic acid conformation: sugar pucker, backbone torsions and
//! a suite sanity check.
//!
//! Nucleotides are recognised by their sugar ring atoms (`C1'`, `C3'`, `C4'`,
//! `O4'`) and measured from atom names alone, so no component dictionary is
//! needed. For each nucleotide the five ring torsions give the pseudorotation
//! phase and amplitude, and the six backbone torsions alpha to zeta plus chi
//! are measured across phosphodiester links.
//!
//! The suite check is a sanity filter, not a named suite clustering. A suite
//! is the backbone stretch between two sugars; for the suite ending at a
//! residue it flags a delta that disagrees with the sugar pucker (on either
//! sugar) and a gamma outside the three staggered windows. Every window comes
//! from the [`NucleicPolicy`] table the reference set resolves, and the result
//! records its identity.
//!
//! Blocks of residues run in parallel and are joined in order, so the arrays
//! are identical at any worker count. Residues the check does not apply to
//! carry `NaN` angles and the [`NOT_EVALUATED`] class code.
//!
//! Complexity: O(residues), constant work per nucleotide (about twenty named
//! atom lookups and a few torsions). Parallel blocks of 256 residues, joined
//! in order.

use molframe_analysis::sugar_pucker;
use molframe_chem::TableIdentity;
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

use super::nucleic_policy::{NucleicPolicy, within};
use super::protein_conformation::{NOT_EVALUATED, to_f32, torsion};
use super::{BuiltinError, Layout, ReferenceSet};

/// Residues per parallel block.
const BLOCK: usize = 256;

/// Pucker class code: C3'-endo.
pub const PUCKER_C3_ENDO: u8 = 0;
/// Pucker class code: C2'-endo.
pub const PUCKER_C2_ENDO: u8 = 1;
/// Pucker class code: neither window.
pub const PUCKER_OTHER: u8 = 2;
/// Suite flag: delta of this sugar disagrees with its pucker.
pub const SUITE_DELTA: u8 = 1;
/// Suite flag: gamma is outside every staggered window.
pub const SUITE_GAMMA: u8 = 2;
/// Suite flag: delta of the preceding sugar disagrees with its pucker.
pub const SUITE_PREVIOUS_DELTA: u8 = 4;
/// Torsions stored per residue: alpha, beta, gamma, delta, epsilon, zeta, chi.
pub const NUCLEIC_TORSIONS: usize = 7;

/// Per-residue nucleic conformation, one entry per structure residue.
#[derive(Clone, Debug, Default)]
pub struct NucleicTable {
    /// Ring torsions nu0 to nu4, row-major, `NaN` when undefined.
    pub nu: Vec<f32>,
    /// Pseudorotation phase in degrees over `[0, 360)`, `NaN` when undefined.
    pub pucker_phase: Vec<f32>,
    /// Pucker amplitude in degrees, `NaN` when undefined.
    pub pucker_amplitude: Vec<f32>,
    /// [`PUCKER_C3_ENDO`], [`PUCKER_C2_ENDO`], [`PUCKER_OTHER`] or
    /// [`NOT_EVALUATED`].
    pub pucker_class: Vec<u8>,
    /// Alpha, beta, gamma, delta, epsilon, zeta, chi per residue
    /// ([`NUCLEIC_TORSIONS`] values), row-major, `NaN` when undefined.
    pub torsions: Vec<f32>,
    /// Bit set of `SUITE_*` flags, or [`NOT_EVALUATED`] without a suite.
    pub suite_flags: Vec<u8>,
    /// 0 sane, 1 flagged, or [`NOT_EVALUATED`] without a suite.
    pub suite_class: Vec<u8>,
    /// The policy the classes were decided by.
    pub tables: Vec<TableIdentity>,
}

impl NucleicTable {
    /// Number of residues covered.
    #[must_use]
    pub fn len(&self) -> usize {
        self.pucker_class.len()
    }

    /// Whether the table is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pucker_class.is_empty()
    }

    fn append(&mut self, other: Self) {
        self.nu.extend(other.nu);
        self.pucker_phase.extend(other.pucker_phase);
        self.pucker_amplitude.extend(other.pucker_amplitude);
        self.pucker_class.extend(other.pucker_class);
        self.torsions.extend(other.torsions);
        self.suite_flags.extend(other.suite_flags);
        self.suite_class.extend(other.suite_class);
    }
}

struct Row {
    nu: [f32; 5],
    phase: f32,
    amplitude: f32,
    class: u8,
    torsions: [f32; NUCLEIC_TORSIONS],
    flags: u8,
    suite: u8,
}

const RING: [[&str; 4]; 5] = [
    ["C4'", "O4'", "C1'", "C2'"],
    ["O4'", "C1'", "C2'", "C3'"],
    ["C1'", "C2'", "C3'", "C4'"],
    ["C2'", "C3'", "C4'", "O4'"],
    ["C3'", "C4'", "O4'", "C1'"],
];

fn ring(layout: &Layout<'_>, residue: usize) -> [Option<f64>; 5] {
    RING.map(|names| {
        torsion(
            layout,
            [
                (residue, names[0]),
                (residue, names[1]),
                (residue, names[2]),
                (residue, names[3]),
            ],
        )
    })
}

fn glycosidic(layout: &Layout<'_>, residue: usize) -> Option<f64> {
    let purine = layout.position(residue, "N9").is_some();
    let (base, ring_atom) = if purine { ("N9", "C4") } else { ("N1", "C2") };
    torsion(
        layout,
        [
            (residue, "O4'"),
            (residue, "C1'"),
            (residue, base),
            (residue, ring_atom),
        ],
    )
}

fn backbone(layout: &Layout<'_>, residue: usize) -> [Option<f64>; 6] {
    let linked_before = layout.phosphate_linked(residue);
    let linked_after = layout.phosphate_linked(residue + 1);
    let next = residue + 1;
    let alpha = match residue.checked_sub(1) {
        Some(previous) if linked_before => torsion(
            layout,
            [
                (previous, "O3'"),
                (residue, "P"),
                (residue, "O5'"),
                (residue, "C5'"),
            ],
        ),
        _ => None,
    };
    let own = |names: [&'static str; 4]| {
        torsion(
            layout,
            [
                (residue, names[0]),
                (residue, names[1]),
                (residue, names[2]),
                (residue, names[3]),
            ],
        )
    };
    let epsilon = linked_after
        .then(|| {
            torsion(
                layout,
                [
                    (residue, "C4'"),
                    (residue, "C3'"),
                    (residue, "O3'"),
                    (next, "P"),
                ],
            )
        })
        .flatten();
    let zeta = linked_after
        .then(|| {
            torsion(
                layout,
                [
                    (residue, "C3'"),
                    (residue, "O3'"),
                    (next, "P"),
                    (next, "O5'"),
                ],
            )
        })
        .flatten();
    [
        alpha,
        own(["P", "O5'", "C5'", "C4'"]),
        own(["O5'", "C5'", "C4'", "C3'"]),
        own(["C5'", "C4'", "C3'", "O3'"]),
        epsilon,
        zeta,
    ]
}

fn sanity(
    policy: &NucleicPolicy,
    gamma: Option<f64>,
    delta: Option<f64>,
    class: u8,
    previous: Option<(u8, f64)>,
) -> (u8, u8) {
    let (Some(gamma), Some(delta), Some((before_class, before_delta))) = (gamma, delta, previous)
    else {
        return (NOT_EVALUATED, NOT_EVALUATED);
    };
    let mut flags = 0u8;
    if !policy.delta_agrees(class, delta) {
        flags |= SUITE_DELTA;
    }
    if !policy.gamma.iter().any(|window| within(gamma, *window)) {
        flags |= SUITE_GAMMA;
    }
    if !policy.delta_agrees(before_class, before_delta) {
        flags |= SUITE_PREVIOUS_DELTA;
    }
    (flags, u8::from(flags != 0))
}

/// Ring torsions nu0 to nu4, each present when its atoms are.
type RingTorsions = [Option<f64>; 5];
/// Phase, amplitude and class of a sugar ring.
type Pucker = (f64, f64, u8);

fn pucker_of(
    policy: &NucleicPolicy,
    layout: &Layout<'_>,
    residue: usize,
) -> (RingTorsions, Option<Pucker>) {
    let nu = ring(layout, residue);
    let complete = match nu {
        [Some(n0), Some(n1), Some(n2), Some(n3), Some(n4)] => Some([n0, n1, n2, n3, n4]),
        _ => None,
    };
    let pucker = complete.map(|angles| {
        let found = sugar_pucker(angles);
        (
            found.phase_degrees,
            found.amplitude,
            policy.pucker_class(found.phase_degrees),
        )
    });
    (nu, pucker)
}

fn previous_sugar(
    policy: &NucleicPolicy,
    layout: &Layout<'_>,
    residue: usize,
) -> Option<(u8, f64)> {
    let before = residue.checked_sub(1)?;
    if !layout.phosphate_linked(residue) || !layout.is_nucleotide(before) {
        return None;
    }
    let (_, pucker) = pucker_of(policy, layout, before);
    let (_, _, class) = pucker?;
    let delta = torsion(
        layout,
        [
            (before, "C5'"),
            (before, "C4'"),
            (before, "C3'"),
            (before, "O3'"),
        ],
    )?;
    Some((class, delta))
}

fn row(policy: &NucleicPolicy, layout: &Layout<'_>, residue: usize) -> Row {
    let mut out = Row {
        nu: [f32::NAN; 5],
        phase: f32::NAN,
        amplitude: f32::NAN,
        class: NOT_EVALUATED,
        torsions: [f32::NAN; NUCLEIC_TORSIONS],
        flags: NOT_EVALUATED,
        suite: NOT_EVALUATED,
    };
    if !layout.is_nucleotide(residue) {
        return out;
    }
    let (nu, pucker) = pucker_of(policy, layout, residue);
    for (slot, value) in out.nu.iter_mut().zip(nu) {
        *slot = to_f32(value);
    }
    if let Some((phase, amplitude, class)) = pucker {
        out.phase = to_f32(Some(phase));
        out.amplitude = to_f32(Some(amplitude));
        out.class = class;
    }
    let angles = backbone(layout, residue);
    for (slot, value) in out.torsions.iter_mut().zip(angles) {
        *slot = to_f32(value);
    }
    out.torsions[6] = to_f32(glycosidic(layout, residue));
    let previous = previous_sugar(policy, layout, residue);
    if out.class != NOT_EVALUATED {
        let (flags, suite) = sanity(policy, angles[2], angles[3], out.class, previous);
        out.flags = flags;
        out.suite = suite;
    }
    out
}

/// Computes the per-residue nucleic conformation of a laid-out structure.
///
/// Runs in `O(residues)` time split across the context's workers. Every window
/// comes from the policy table of `reference`, whose identity the result records.
///
/// # Errors
///
/// Returns [`BuiltinError`] if no policy table is selected or it is bad, or a
/// worker panics.
pub fn nucleic_table(
    layout: &Layout<'_>,
    reference: &ReferenceSet,
    context: &ExecutionContext,
) -> Result<NucleicTable, BuiltinError> {
    let resolved = reference.nucleic()?;
    let policy = &resolved.table;
    let plan = BlockPlan::new(layout.len(), BLOCK);
    let blocks = map_blocks_in(plan, context, |_, range| {
        let mut part = NucleicTable::default();
        for residue in range {
            let found = row(policy, layout, residue);
            part.nu.extend_from_slice(&found.nu);
            part.pucker_phase.push(found.phase);
            part.pucker_amplitude.push(found.amplitude);
            part.pucker_class.push(found.class);
            part.torsions.extend_from_slice(&found.torsions);
            part.suite_flags.push(found.flags);
            part.suite_class.push(found.suite);
        }
        part
    })?;
    let mut table = NucleicTable {
        tables: vec![resolved.identity.clone()],
        ..NucleicTable::default()
    };
    for block in blocks {
        table.append(block);
    }
    Ok(table)
}

#[cfg(test)]
#[path = "nucleic_scoring_tests.rs"]
mod tests;
