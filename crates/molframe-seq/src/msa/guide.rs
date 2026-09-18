//! Linear-storage single-linkage guide construction.

use super::MsaError;
use super::profile::GAP;

pub(super) fn single_linkage_guide(
    sequences: &[&[u8]],
    memory_limit_bytes: usize,
) -> Result<Vec<(usize, usize)>, MsaError> {
    if sequences.len() < 2 {
        return Ok(Vec::new());
    }
    let size = sequences.len();
    let max_length = match sequences
        .iter()
        .map(|sequence| ungapped_len(sequence))
        .max()
    {
        Some(length) => length,
        None => 0,
    };
    let (symbol_slots, alphabet_size) = symbol_slots(sequences);
    let required = guide_workspace_bytes(size, max_length, alphabet_size)?;
    check_memory(required, memory_limit_bytes)?;
    let mut lengths = Vec::new();
    lengths
        .try_reserve_exact(size)
        .map_err(|_| MsaError::MemoryLimit {
            required,
            limit: memory_limit_bytes,
        })?;
    lengths.extend(sequences.iter().map(|sequence| ungapped_len(sequence)));
    let mut lcs = LcsWorkspace::new(
        max_length,
        &symbol_slots,
        alphabet_size,
        required,
        memory_limit_bytes,
    )?;
    let mut selected = vec![false; size];
    let mut distances = vec![f64::INFINITY; size];
    let mut nearest = vec![0_usize; size];
    let mut edges = Vec::with_capacity(size - 1);
    selected[0] = true;
    for candidate in 1..size {
        distances[candidate] = lcs.distance_with_lengths(
            sequences[0],
            lengths[0],
            sequences[candidate],
            lengths[candidate],
        )?;
    }
    for _ in 1..size {
        let next = (0..size)
            .filter(|candidate| !selected[*candidate])
            .min_by(|left, right| {
                distances[*left]
                    .total_cmp(&distances[*right])
                    .then_with(|| left.cmp(right))
            })
            .ok_or(MsaError::DimensionOverflow)?;
        edges.push(MstEdge {
            distance: distances[next],
            first: nearest[next],
            second: next,
        });
        selected[next] = true;
        for candidate in 0..size {
            if selected[candidate] {
                continue;
            }
            let distance = lcs.distance_with_lengths(
                sequences[next],
                lengths[next],
                sequences[candidate],
                lengths[candidate],
            )?;
            let ordering = distance.total_cmp(&distances[candidate]);
            if ordering.is_lt() || (ordering.is_eq() && next < nearest[candidate]) {
                distances[candidate] = distance;
                nearest[candidate] = next;
            }
        }
    }
    edges.sort_by(|left, right| {
        left.distance
            .total_cmp(&right.distance)
            .then_with(|| left.ordered().cmp(&right.ordered()))
    });
    Ok(merges_from_mst(&edges, size))
}

struct LcsWorkspace {
    word_capacity: usize,
    masks: Vec<u64>,
    state: Vec<u64>,
    touched: Vec<u8>,
    seen: [bool; 256],
    symbol_slots: [u16; 256],
}

impl LcsWorkspace {
    fn new(
        max_length: usize,
        symbol_slots: &[u16; 256],
        alphabet_size: usize,
        required: usize,
        limit: usize,
    ) -> Result<Self, MsaError> {
        let word_capacity = max_length
            .checked_add(63)
            .map(|length| length / 64)
            .ok_or(MsaError::DimensionOverflow)?;
        let mask_words = word_capacity
            .checked_mul(alphabet_size)
            .ok_or(MsaError::DimensionOverflow)?;
        let mut masks = Vec::new();
        masks
            .try_reserve_exact(mask_words)
            .map_err(|_| MsaError::MemoryLimit { required, limit })?;
        masks.resize(mask_words, 0);
        let mut state = Vec::new();
        state
            .try_reserve_exact(word_capacity)
            .map_err(|_| MsaError::MemoryLimit { required, limit })?;
        state.resize(word_capacity, 0);
        let mut touched = Vec::new();
        touched
            .try_reserve_exact(alphabet_size)
            .map_err(|_| MsaError::MemoryLimit { required, limit })?;
        Ok(Self {
            word_capacity,
            masks,
            state,
            touched,
            seen: [false; 256],
            symbol_slots: *symbol_slots,
        })
    }

    #[cfg(test)]
    fn distance(&mut self, first: &[u8], second: &[u8]) -> Result<f64, MsaError> {
        self.distance_with_lengths(first, ungapped_len(first), second, ungapped_len(second))
    }

    fn distance_with_lengths(
        &mut self,
        first: &[u8],
        first_length: usize,
        second: &[u8],
        second_length: usize,
    ) -> Result<f64, MsaError> {
        let (text, text_length, pattern, pattern_length) = if first_length >= second_length {
            (first, first_length, second, second_length)
        } else {
            (second, second_length, first, first_length)
        };
        let normalizer = pattern_length.max(text_length);
        if normalizer == 0 {
            return Ok(0.0);
        }
        if pattern_length == 0 {
            return Ok(1.0);
        }
        let words = pattern_length.div_ceil(64);
        self.clear_masks();
        let mut position = 0_usize;
        for &symbol in pattern {
            if symbol == GAP {
                continue;
            }
            let symbol_index = usize::from(symbol);
            if !self.seen[symbol_index] {
                self.seen[symbol_index] = true;
                self.touched.push(symbol);
            }
            let slot = usize::from(self.symbol_slots[symbol_index]);
            let word = position / 64;
            let bit = position % 64;
            self.masks[slot * self.word_capacity + word] |= 1_u64 << bit;
            position += 1;
        }
        self.state[..words].fill(0);
        for &symbol in text {
            if symbol == GAP {
                continue;
            }
            let slot = self.symbol_slots[usize::from(symbol)];
            if slot == u16::MAX {
                continue;
            }
            let mask_start = usize::from(slot) * self.word_capacity;
            let mut shift_carry = 1_u64;
            let mut borrow = false;
            for word in 0..words {
                let previous = self.state[word];
                let matches = self.masks[mask_start + word];
                let union = matches | previous;
                let shifted = (previous << 1) | shift_carry;
                shift_carry = previous >> 63;
                let (partial, first_borrow) = union.overflowing_sub(shifted);
                let (difference, second_borrow) = partial.overflowing_sub(u64::from(borrow));
                borrow = first_borrow || second_borrow;
                self.state[word] = union & !difference;
            }
        }
        let common = self.state[..words]
            .iter()
            .try_fold(0_u32, |total, word| total.checked_add(word.count_ones()))
            .ok_or(MsaError::DimensionOverflow)?;
        let normalizer = u32::try_from(normalizer).map_err(|_| MsaError::DimensionOverflow)?;
        Ok(1.0 - f64::from(common) / f64::from(normalizer))
    }

    fn clear_masks(&mut self) {
        for symbol in self.touched.drain(..) {
            self.seen[usize::from(symbol)] = false;
            let slot = usize::from(self.symbol_slots[usize::from(symbol)]);
            let start = slot * self.word_capacity;
            self.masks[start..start + self.word_capacity].fill(0);
        }
    }
}

fn ungapped_len(sequence: &[u8]) -> usize {
    sequence.iter().filter(|symbol| **symbol != GAP).count()
}

fn guide_workspace_bytes(
    size: usize,
    max_length: usize,
    alphabet_size: usize,
) -> Result<usize, MsaError> {
    let words = max_length
        .checked_add(63)
        .map(|length| length / 64)
        .ok_or(MsaError::DimensionOverflow)?;
    let lcs = words
        .checked_mul(
            alphabet_size
                .checked_add(1)
                .ok_or(MsaError::DimensionOverflow)?,
        )
        .and_then(|words| words.checked_mul(size_of::<u64>()))
        .ok_or(MsaError::DimensionOverflow)?;
    let per_site =
        size_of::<bool>() + size_of::<f64>() + 6 * size_of::<usize>() + size_of::<MstEdge>();
    size.checked_mul(per_site)
        .and_then(|guide| guide.checked_add(lcs))
        .and_then(|bytes| bytes.checked_add(alphabet_size + 3 * 256))
        .ok_or(MsaError::DimensionOverflow)
}

fn symbol_slots(sequences: &[&[u8]]) -> ([u16; 256], usize) {
    let mut slots = [u16::MAX; 256];
    let mut next = 0u16;
    for symbol in sequences
        .iter()
        .flat_map(|sequence| sequence.iter().copied())
        .filter(|symbol| *symbol != GAP)
    {
        let slot = &mut slots[usize::from(symbol)];
        if *slot == u16::MAX {
            *slot = next;
            next += 1;
        }
    }
    (slots, usize::from(next))
}

fn check_memory(required: usize, limit: usize) -> Result<(), MsaError> {
    if required <= limit {
        Ok(())
    } else {
        Err(MsaError::MemoryLimit { required, limit })
    }
}

#[derive(Clone, Copy)]
struct MstEdge {
    distance: f64,
    first: usize,
    second: usize,
}

impl MstEdge {
    fn ordered(self) -> (usize, usize) {
        if self.first < self.second {
            (self.first, self.second)
        } else {
            (self.second, self.first)
        }
    }
}

fn merges_from_mst(edges: &[MstEdge], size: usize) -> Vec<(usize, usize)> {
    let mut parent = (0..size).collect::<Vec<_>>();
    let mut cluster_slot = (0..size).collect::<Vec<_>>();
    let mut merges = Vec::with_capacity(edges.len());
    for edge in edges {
        let first_root = find_root(&mut parent, edge.first);
        let second_root = find_root(&mut parent, edge.second);
        if first_root == second_root {
            continue;
        }
        let (left, right) = if cluster_slot[first_root] < cluster_slot[second_root] {
            (cluster_slot[first_root], cluster_slot[second_root])
        } else {
            (cluster_slot[second_root], cluster_slot[first_root])
        };
        merges.push((left, right));
        parent[second_root] = first_root;
        cluster_slot[first_root] = left;
    }
    merges
}

fn find_root(parent: &mut [usize], mut node: usize) -> usize {
    while parent[node] != node {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    node
}

#[cfg(test)]
#[path = "guide_tests.rs"]
mod tests;
