//! Residue numbers written as ranges, and as a `resid` query over them.

use std::collections::BTreeMap;

/// Consecutive runs `(first, last)` of the sorted, distinct `numbers`.
fn runs(numbers: &[i32]) -> Vec<(i32, i32)> {
    let mut out: Vec<(i32, i32)> = Vec::new();
    for &number in numbers {
        match out.last_mut() {
            Some((_, end)) if end.checked_add(1) == Some(number) => *end = number,
            _ => out.push((number, number)),
        }
    }
    out
}

fn sorted_distinct(numbers: &[i32]) -> Vec<i32> {
    let mut sorted = numbers.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    sorted
}

/// `[20, 21, 22, 25]` as `20-22,25`. The input need not be sorted and may repeat.
#[must_use]
pub fn compact_ranges(numbers: &[i32]) -> String {
    runs(&sorted_distinct(numbers))
        .into_iter()
        .map(|(start, end)| {
            if start == end {
                start.to_string()
            } else {
                format!("{start}-{end}")
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// A `resid` query for residues given as (chain, number) pairs: consecutive numbers share a
/// range (`resid 3:7`), and each chain is its own scope (an empty chain name has none). Returns
/// the query and the number of distinct residues, or `None` when there are none.
#[must_use]
pub fn resid_query(residues: &[(String, i32)]) -> Option<(String, usize)> {
    let mut by_chain: BTreeMap<&str, Vec<i32>> = BTreeMap::new();
    for (chain, number) in residues {
        by_chain.entry(chain).or_default().push(*number);
    }
    let mut parts = Vec::new();
    let mut count = 0;
    for (chain, numbers) in &by_chain {
        let numbers = sorted_distinct(numbers);
        count += numbers.len();
        let joined = runs(&numbers)
            .into_iter()
            .map(|(start, end)| {
                if start == end {
                    format!("resid {start}")
                } else {
                    format!("resid {start}:{end}")
                }
            })
            .collect::<Vec<_>>()
            .join(" or ");
        parts.push(if chain.is_empty() {
            format!("({joined})")
        } else {
            format!("(chain {chain} and ({joined}))")
        });
    }
    (count > 0).then(|| (parts.join(" or "), count))
}

#[cfg(test)]
#[path = "residue_ranges_tests.rs"]
mod tests;
