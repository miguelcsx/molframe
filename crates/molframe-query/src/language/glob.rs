//! Deterministic glob matching without a regular-expression engine.

#[derive(Clone, PartialEq, Eq, Debug)]
enum Part {
    Literal(char),
    One,
    Any,
    Class {
        negated: bool,
        ranges: Vec<(u32, u32)>,
    },
}

/// A compiled identifier glob.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) struct Glob {
    parts: Vec<Part>,
}

impl Glob {
    /// Compiles `pattern` into deterministic glob components.
    ///
    /// Compilation is `O(P log R)` in the worst case due to normalization of
    /// character-class ranges. Unlike the original implementation, Unicode
    /// ranges such as `[a-z]` are stored as intervals rather than expanded into
    /// one `char` per member.
    pub(crate) fn new(pattern: &str) -> Self {
        let mut parts = Vec::new();
        let mut position = 0usize;

        while let Some((character, next)) = character_at(pattern, position) {
            match character {
                '*' => {
                    if !matches!(parts.last(), Some(Part::Any)) {
                        parts.push(Part::Any);
                    }

                    position = next;
                }
                '?' => {
                    parts.push(Part::One);
                    position = next;
                }
                '[' => {
                    if let Some((part, next)) = class(pattern, position) {
                        parts.push(part);
                        position = next;
                    } else {
                        parts.push(Part::Literal('['));
                        position = next;
                    }
                }
                '\\' => {
                    if let Some((escaped, escaped_next)) = character_at(pattern, next) {
                        parts.push(Part::Literal(escaped));
                        position = escaped_next;
                    } else {
                        parts.push(Part::Literal('\\'));
                        position = next;
                    }
                }
                literal => {
                    parts.push(Part::Literal(literal));
                    position = next;
                }
            }
        }

        Self { parts }
    }

    /// Returns whether `text` satisfies this glob.
    ///
    /// Matching uses the classic greedy wildcard algorithm with backtracking to
    /// the latest `*`. It performs no heap allocation. Typical runtime is
    /// linear in pattern plus text length; worst-case runtime remains
    /// `O(P * T)` for adversarial wildcard patterns, while auxiliary space is
    /// reduced from `O(T)` to `O(1)`.
    pub(crate) fn matches(&self, text: &str) -> bool {
        let mut part_position = 0usize;
        let mut text_position = 0usize;
        let mut star_part = None;
        let mut star_text = 0usize;

        while text_position < text.len() {
            let Some((character, next_text)) = character_at(text, text_position) else {
                return false;
            };

            match self.parts.get(part_position) {
                Some(Part::Any) => {
                    star_part = Some(part_position);
                    star_text = text_position;
                    part_position += 1;
                }
                Some(part) if part_matches(part, character) => {
                    part_position += 1;
                    text_position = next_text;
                }
                _ => {
                    let Some(star) = star_part else {
                        return false;
                    };

                    let Some((_, next_star_text)) = character_at(text, star_text) else {
                        return false;
                    };

                    star_text = next_star_text;
                    text_position = next_star_text;
                    part_position = star + 1;
                }
            }
        }

        while matches!(self.parts.get(part_position), Some(Part::Any)) {
            part_position += 1;
        }

        part_position == self.parts.len()
    }
}

/// Parses a character class beginning at byte offset `start`.
///
/// Members and ranges are normalized into sorted non-overlapping Unicode scalar
/// intervals. Returns `None` for an unterminated class.
fn class(pattern: &str, start: usize) -> Option<(Part, usize)> {
    let mut position = start.checked_add(1)?;
    let mut negated = false;

    if character_at(pattern, position).is_some_and(|(character, _)| character == '!') {
        let (_, next) = character_at(pattern, position)?;
        negated = true;
        position = next;
    }

    let mut ranges = Vec::new();

    while let Some((character, next)) = character_at(pattern, position) {
        if character == ']' && !ranges.is_empty() {
            normalize_ranges(&mut ranges);

            return Some((Part::Class { negated, ranges }, next));
        }

        if let Some(('-', after_dash)) = character_at(pattern, next)
            && let Some((end, after_end)) = character_at(pattern, after_dash)
            && end != ']'
        {
            let start_code = u32::from(character);
            let end_code = u32::from(end);

            ranges.push((start_code.min(end_code), start_code.max(end_code)));

            position = after_end;
            continue;
        }

        let code = u32::from(character);
        ranges.push((code, code));
        position = next;
    }

    None
}

/// Sorts and merges overlapping or adjacent character-class intervals.
///
/// Runtime is `O(R log R)` for `R` input ranges and additional space is `O(1)`
/// beyond the existing vector.
fn normalize_ranges(ranges: &mut Vec<(u32, u32)>) {
    if ranges.len() < 2 {
        return;
    }

    ranges.sort_unstable();

    let mut write = 0usize;

    for read in 1..ranges.len() {
        let (start, end) = ranges[read];
        let current_end = ranges[write].1;

        let adjacent_or_overlapping = match current_end.checked_add(1) {
            Some(boundary) => start <= boundary,
            None => true,
        };
        if adjacent_or_overlapping {
            if end > current_end {
                ranges[write].1 = end;
            }
        } else {
            write += 1;
            ranges[write] = (start, end);
        }
    }

    ranges.truncate(write + 1);
}

/// Returns the Unicode scalar beginning at byte offset `position` and its end.
///
/// Invalid offsets return `None`; no slicing panic is possible.
#[inline]
fn character_at(text: &str, position: usize) -> Option<(char, usize)> {
    let tail = text.get(position..)?;
    let character = tail.chars().next()?;

    Some((character, position.checked_add(character.len_utf8())?))
}

/// Tests one compiled glob part against a Unicode scalar.
///
/// Character classes use binary search over normalized intervals, so class
/// matching costs `O(log R)` rather than a linear scan over expanded members.
#[inline]
fn part_matches(part: &Part, character: char) -> bool {
    match part {
        Part::Literal(expected) => *expected == character,
        Part::One | Part::Any => true,
        Part::Class { negated, ranges } => {
            let code = u32::from(character);
            let position = ranges.partition_point(|(start, _)| *start <= code);

            let contained = position
                .checked_sub(1)
                .and_then(|index| ranges.get(index))
                .is_some_and(|(_, end)| code <= *end);

            contained != *negated
        }
    }
}

#[cfg(test)]
#[path = "glob_tests.rs"]
mod tests;
