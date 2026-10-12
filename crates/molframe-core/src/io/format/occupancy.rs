//! The read-time occupancy policy.

use super::options::OccupancyPolicy;
use crate::chunk::ChunkBuilder;
use crate::diagnostic::{Code, Diagnostic, Diagnostics};
use crate::structure::{CoordinateStore, Structure};
use std::sync::Arc;

#[derive(Default)]
struct Scan {
    count: u64,
    first: Option<(u32, f32)>,
}

impl Scan {
    fn add(&mut self, other: Self) {
        self.count += other.count;
        if self.first.is_none() {
            self.first = other.first;
        }
    }
}

fn scan(structure: &Structure) -> Scan {
    let mut found = Scan::default();
    if let Some(models) = structure.ragged_models() {
        for model in models {
            found.add(scan(model));
        }
        return found;
    }
    for chunk in structure.data().chunks.iter() {
        for local in 0..chunk.len() {
            let Some((value, presence)) = chunk.occupancy(local) else {
                continue;
            };
            if !presence.is_present() || (0.0..=1.0).contains(&value) {
                continue;
            }
            found.count += 1;
            if found.first.is_none() {
                found.first = Some((chunk.atoms().start.saturating_add(local), value));
            }
        }
    }
    found
}

fn clamped(structure: &Structure) -> Option<Structure> {
    let data = structure.data();
    let mut next = data.clone();
    if let Some(models) = structure.ragged_models() {
        let models = models
            .iter()
            .map(|m| clamped(m))
            .collect::<Option<Vec<_>>>()?;
        next.coords = CoordinateStore::Ragged { models };
        return Some(Structure::new(next));
    }
    let mut builder = ChunkBuilder::new();
    for chunk in data.chunks.iter() {
        builder.start_model(chunk.model());
        for atom in chunk.atoms() {
            let local = atom.checked_sub(chunk.atoms().start)?;
            let mut record = chunk.record(local, &data.topology.residues, Some([0.0; 3]))?;
            if record.occupancy.1.is_present() {
                record.occupancy.0 = record.occupancy.0.clamp(0.0, 1.0);
            }
            builder.push(record);
        }
    }
    next.chunks = Arc::new(builder.finish().0);
    Some(Structure::new(next))
}

/// Applies the policy; returns the structure and whether the read must fail.
pub(super) fn apply(
    structure: Structure,
    policy: OccupancyPolicy,
    findings: &mut Diagnostics,
) -> (Structure, bool) {
    let found = scan(&structure);
    let Some((atom, value)) = found.first else {
        return (structure, false);
    };
    let context = |code| {
        Diagnostic::new(code)
            .with_context("atoms", found.count.to_string())
            .with_context("first_atom", atom.to_string())
            .with_context("first_occupancy", value.to_string())
    };
    match policy {
        OccupancyPolicy::Keep => {
            findings.push(context(Code::E3009));
            (structure, false)
        }
        OccupancyPolicy::Reject => {
            findings.push(context(Code::E3009));
            (structure, true)
        }
        OccupancyPolicy::Clamp => match clamped(&structure) {
            Some(next) => {
                findings.push(context(Code::W3014));
                (next, false)
            }
            None => {
                findings.push(context(Code::E3009));
                (structure, true)
            }
        },
    }
}

#[cfg(test)]
#[path = "occupancy_tests.rs"]
mod tests;
