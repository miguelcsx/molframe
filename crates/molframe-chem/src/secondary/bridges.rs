//! Beta bridges, ladders and their continuous bulge-spanning strand intervals.

use super::geometry::continuous;
use super::{Bonds, DsspBackbone as Backbone};
use molframe_core::ExecutionContext;
use molframe_core::SecondaryStructure as Ss;
use molframe_core::hashing::IdentityHashMap;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Orientation {
    Parallel,
    Antiparallel,
}

#[derive(Clone, Copy)]
struct Ladder {
    first: usize,
    last: usize,
    low: usize,
    high: usize,
    bridges: usize,
    orientation: Orientation,
}

const PAIRS_PER_BLOCK: usize = 65_536;

/// The pairs that form a bridge, with their orientation, in pair order.
///
/// The test reads only the backbones and the hydrogen bonds, so it runs in
/// blocks; the ladder fold that follows is order dependent and stays serial.
fn bridge_candidates(
    backbones: &[Backbone],
    bonds: &Bonds,
    pairs: &[(usize, usize)],
    context: &ExecutionContext,
) -> Vec<(usize, usize, Orientation)> {
    let has = |a, b| bonds.contains(&(a, b));
    let scan = |range: std::ops::Range<usize>| {
        let mut found = Vec::new();
        for &(i, j) in &pairs[range] {
            if i == 0
                || j == 0
                || !continuous(backbones, i - 1, i + 1)
                || !continuous(backbones, j - 1, j + 1)
            {
                continue;
            }
            let orientation =
                if (has(i - 1, j) && has(j, i + 1)) || (has(j - 1, i) && has(i, j + 1)) {
                    Orientation::Parallel
                } else if (has(i, j) && has(j, i)) || (has(i - 1, j + 1) && has(j - 1, i + 1)) {
                    Orientation::Antiparallel
                } else {
                    continue;
                };
            found.push((i, j, orientation));
        }
        found
    };
    let plan = BlockPlan::new(pairs.len(), PAIRS_PER_BLOCK);
    let blocks = match map_blocks_in(plan, context, |_, range| scan(range)) {
        Ok(blocks) => blocks,
        Err(_) => vec![scan(0..pairs.len())],
    };
    blocks.into_iter().flatten().collect()
}

pub(super) fn assign(
    backbones: &[Backbone],
    bonds: &Bonds,
    pairs: &[(usize, usize)],
    states: &mut [Ss],
    context: &ExecutionContext,
) {
    let mut ladders: Vec<Ladder> = Vec::new();
    let mut ends = IdentityHashMap::default();
    for (i, j, orientation) in bridge_candidates(backbones, bonds, pairs, context) {
        let prior_j = match orientation {
            Orientation::Parallel => j.checked_sub(1),
            Orientation::Antiparallel => j.checked_add(1),
        };
        let previous = prior_j.and_then(|p| ends.remove(&(i - 1, p, orientation)));
        let index = if let Some(index) = previous {
            let ladder: &mut Ladder = &mut ladders[index];
            ladder.last = i;
            ladder.low = ladder.low.min(j);
            ladder.high = ladder.high.max(j);
            ladder.bridges += 1;
            index
        } else {
            let index = ladders.len();
            ladders.push(Ladder {
                first: i,
                last: i,
                low: j,
                high: j,
                bridges: 1,
                orientation,
            });
            index
        };
        ends.insert((i, j, orientation), index);
    }
    // Ladders are already ordered by the first arm. Only a five-residue local
    // neighbourhood can join across a bulge; no all-pairs residue scan is needed.
    for i in 0..ladders.len() {
        if ladders[i].bridges == 0 {
            continue;
        }
        for j in i + 1..ladders.len() {
            let (a, b) = (ladders[i], ladders[j]);
            if b.first >= a.last.saturating_add(6) {
                break;
            }
            if b.bridges == 0 || b.orientation != a.orientation || b.first <= a.last {
                continue;
            }
            let gap = match a.orientation {
                Orientation::Parallel => b.low.checked_sub(a.high),
                Orientation::Antiparallel => a.low.checked_sub(b.high),
            };
            let Some(gap) = gap else {
                continue;
            };
            if !((gap < 6 && b.first - a.last < 3) || gap < 3)
                || !continuous(backbones, a.first, b.last)
                || !continuous(backbones, a.low.min(b.low), a.high.max(b.high))
            {
                continue;
            }
            ladders[i].last = b.last;
            ladders[i].low = a.low.min(b.low);
            ladders[i].high = a.high.max(b.high);
            ladders[i].bridges += b.bridges;
            ladders[j].bridges = 0;
        }
    }
    for ladder in ladders {
        if ladder.bridges == 0 {
            continue;
        }
        let kind = if ladder.bridges > 1 {
            Ss::Strand
        } else {
            Ss::BetaBridge
        };
        for residue in (ladder.first..=ladder.last).chain(ladder.low..=ladder.high) {
            if states[residue] != Ss::Strand {
                states[residue] = kind;
            }
        }
    }
}
