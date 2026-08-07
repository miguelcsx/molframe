//! Deterministic glob matching without a regular-expression engine.

#[derive(Clone, PartialEq, Eq, Debug)]
enum Part {
    Literal(char),
    One,
    Any,
    Class { negated: bool, members: Vec<char> },
}

/// A compiled identifier glob.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(super) struct Glob {
    parts: Vec<Part>,
}

impl Glob {
    pub(super) fn new(pattern: &str) -> Self {
        let characters: Vec<char> = pattern.chars().collect();
        let mut parts = Vec::new();
        let mut position = 0usize;
        while position < characters.len() {
            match characters[position] {
                '*' => {
                    if !matches!(parts.last(), Some(Part::Any)) {
                        parts.push(Part::Any);
                    }
                    position += 1;
                }
                '?' => {
                    parts.push(Part::One);
                    position += 1;
                }
                '[' => {
                    if let Some((part, next)) = class(&characters, position) {
                        parts.push(part);
                        position = next;
                    } else {
                        parts.push(Part::Literal('['));
                        position += 1;
                    }
                }
                '\\' if position + 1 < characters.len() => {
                    parts.push(Part::Literal(characters[position + 1]));
                    position += 2;
                }
                literal => {
                    parts.push(Part::Literal(literal));
                    position += 1;
                }
            }
        }
        Self { parts }
    }

    pub(super) fn matches(&self, text: &str) -> bool {
        let characters: Vec<char> = text.chars().collect();
        let mut previous = vec![false; characters.len() + 1];
        previous[0] = true;
        for part in &self.parts {
            let mut current = vec![false; characters.len() + 1];
            match part {
                Part::Any => {
                    current[0] = previous[0];
                    for position in 1..=characters.len() {
                        current[position] = previous[position] || current[position - 1];
                    }
                }
                _ => {
                    for position in 1..=characters.len() {
                        current[position] =
                            previous[position - 1] && part_matches(part, characters[position - 1]);
                    }
                }
            }
            previous = current;
        }
        previous.last().copied().is_some_and(|matched| matched)
    }
}

fn class(characters: &[char], start: usize) -> Option<(Part, usize)> {
    let mut position = start + 1;
    let negated = characters
        .get(position)
        .is_some_and(|character| *character == '!');
    if negated {
        position += 1;
    }
    let mut members = Vec::new();
    while let Some(character) = characters.get(position).copied() {
        if character == ']' && !members.is_empty() {
            return Some((Part::Class { negated, members }, position + 1));
        }
        if position + 2 < characters.len() && characters[position + 1] == '-' {
            let end = characters[position + 2];
            if end != ']' {
                let start_code = u32::from(character);
                let end_code = u32::from(end);
                for code in start_code.min(end_code)..=start_code.max(end_code) {
                    if let Some(member) = char::from_u32(code) {
                        members.push(member);
                    }
                }
                position += 3;
                continue;
            }
        }
        members.push(character);
        position += 1;
    }
    None
}

fn part_matches(part: &Part, character: char) -> bool {
    match part {
        Part::Literal(expected) => *expected == character,
        Part::One | Part::Any => true,
        Part::Class { negated, members } => members.contains(&character) != *negated,
    }
}

#[cfg(test)]
#[path = "glob_tests.rs"]
mod tests;
