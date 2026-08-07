//! Python sequence-index conversion at the language boundary.

pub(crate) fn normalise_index(index: isize, count: usize) -> Option<usize> {
    let count = isize::try_from(count).ok()?;
    let position = if index < 0 {
        count.checked_add(index)?
    } else {
        index
    };
    if !(0..count).contains(&position) {
        return None;
    }
    usize::try_from(position).ok()
}

#[cfg(test)]
#[path = "index_tests.rs"]
mod tests;
