use crate::StereoConfiguration;
use crate::smarts::{
    AtomExpression, AtomTest, BondExpression, PatternBond, SignedAtomTest, SmartsError,
    SmartsPattern,
};
use molframe_core::{BondOrder, Element};
use std::collections::BTreeMap;

#[path = "smarts_parse/nesting.rs"]
mod nesting;
use nesting::split_top_level;

pub(crate) fn parse(text: &str) -> Result<SmartsPattern, SmartsError> {
    Parser::new(text).parse()
}

struct Parser<'a> {
    text: &'a str,
    bytes: &'a [u8],
    position: usize,
    atoms: Vec<AtomExpression>,
    bonds: Vec<PatternBond>,
    current: Option<usize>,
    branches: Vec<usize>,
    rings: BTreeMap<u16, (usize, Option<BondExpression>)>,
    pending_bond: Option<BondExpression>,
}

impl<'a> Parser<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            text,
            bytes: text.as_bytes(),
            position: 0,
            atoms: Vec::new(),
            bonds: Vec::new(),
            current: None,
            branches: Vec::new(),
            rings: BTreeMap::new(),
            pending_bond: None,
        }
    }

    fn parse(mut self) -> Result<SmartsPattern, SmartsError> {
        if self.bytes.is_empty() {
            return Err(self.error("empty pattern"));
        }
        while self.position < self.bytes.len() {
            match self.bytes[self.position] {
                b'(' => self.open_branch()?,
                b')' => self.close_branch()?,
                b'.' => self.disconnect()?,
                b'-' | b'=' | b'#' | b':' | b'~' | b'@' | b'/' | b'\\' | b'!' => {
                    self.read_bond()?;
                }
                b'0'..=b'9' | b'%' => self.read_ring()?,
                _ => self.read_atom()?,
            }
        }
        if !self.branches.is_empty() {
            return Err(self.error("unclosed branch"));
        }
        if !self.rings.is_empty() {
            return Err(self.error("unclosed ring bond"));
        }
        if self.pending_bond.is_some() {
            return Err(self.error("bond without a following atom or ring"));
        }
        Ok(SmartsPattern {
            atoms: self.atoms,
            bonds: self.bonds,
        })
    }

    fn open_branch(&mut self) -> Result<(), SmartsError> {
        let current = self
            .current
            .ok_or_else(|| self.error("branch has no parent atom"))?;
        self.branches.push(current);
        self.position += 1;
        Ok(())
    }

    fn close_branch(&mut self) -> Result<(), SmartsError> {
        if self.pending_bond.is_some() {
            return Err(self.error("branch ends after a bond"));
        }
        self.current = Some(
            self.branches
                .pop()
                .ok_or_else(|| self.error("unmatched closing parenthesis"))?,
        );
        self.position += 1;
        Ok(())
    }

    fn disconnect(&mut self) -> Result<(), SmartsError> {
        if self.current.is_none() || self.pending_bond.is_some() {
            return Err(self.error("invalid disconnected-component separator"));
        }
        self.current = None;
        self.position += 1;
        Ok(())
    }

    fn read_atom(&mut self) -> Result<(), SmartsError> {
        let expression = if self.bytes[self.position] == b'[' {
            let start = self.position + 1;
            let end = self.bracket_end(start)?;
            let expression = parse_bracket(&self.text[start..end], start)?;
            self.position = end + 1;
            expression
        } else {
            self.read_bare_atom()?
        };
        let next = self.atoms.len();
        self.atoms.push(expression);
        if let Some(previous) = self.current {
            let expression = match self.pending_bond.take() {
                Some(expression) => expression,
                None => BondExpression::Default,
            };
            self.bonds.push(PatternBond {
                first: previous,
                second: next,
                expression,
            });
        } else if self.pending_bond.is_some() {
            return Err(self.error("bond has no preceding atom"));
        }
        self.current = Some(next);
        Ok(())
    }

    fn bracket_end(&self, start: usize) -> Result<usize, SmartsError> {
        let mut recursive_depth = 0usize;
        for position in start..self.bytes.len() {
            match self.bytes[position] {
                b'(' if position > 0 && self.bytes[position - 1] == b'$' => recursive_depth += 1,
                b'(' if recursive_depth > 0 => recursive_depth += 1,
                b')' if recursive_depth > 0 => recursive_depth -= 1,
                b']' if recursive_depth == 0 => return Ok(position),
                _ => {}
            }
        }
        Err(SmartsError::new(
            if start == 0 { 0 } else { start - 1 },
            "unclosed atom expression",
        ))
    }

    fn read_bare_atom(&mut self) -> Result<AtomExpression, SmartsError> {
        let start = self.position;
        if self.bytes[start] == b'*' {
            self.position += 1;
            return Ok(single(AtomTest::Any));
        }
        let aromatic = matches!(self.bytes[start], b'b' | b'c' | b'n' | b'o' | b'p' | b's');
        let length = if start + 1 < self.bytes.len()
            && self.bytes[start].is_ascii_uppercase()
            && self.bytes[start + 1].is_ascii_lowercase()
        {
            2
        } else {
            1
        };
        let token = &self.text[start..start + length];
        let element = Element::from_symbol(token)
            .ok_or_else(|| SmartsError::new(start, "unknown bare atom symbol"))?;
        self.position += length;
        let mut tests = vec![SignedAtomTest {
            negated: false,
            test: AtomTest::Element(element),
        }];
        if aromatic {
            tests.push(SignedAtomTest {
                negated: false,
                test: AtomTest::Aromatic,
            });
        }
        Ok(AtomExpression {
            alternatives: vec![tests],
        })
    }

    fn read_bond(&mut self) -> Result<(), SmartsError> {
        if self.current.is_none() || self.pending_bond.is_some() {
            return Err(self.error("misplaced or repeated bond expression"));
        }
        let start = self.position;
        let negated = self.bytes[start] == b'!';
        if negated {
            self.position += 1;
        }
        let byte = *self
            .bytes
            .get(self.position)
            .ok_or_else(|| self.error("incomplete negated bond"))?;
        self.position += 1;
        let expression = match (negated, byte) {
            (false, b'-' | b'/' | b'\\') => BondExpression::Order(BondOrder::Single),
            (false, b'=') => BondExpression::Order(BondOrder::Double),
            (false, b'#') => BondExpression::Order(BondOrder::Triple),
            (false, b':') => BondExpression::Order(BondOrder::Aromatic),
            (false, b'~') => BondExpression::Any,
            (false, b'@') => BondExpression::Ring(true),
            (true, b'@') => BondExpression::Ring(false),
            (true, b'-') => BondExpression::NotOrder(BondOrder::Single),
            (true, b'=') => BondExpression::NotOrder(BondOrder::Double),
            (true, b'#') => BondExpression::NotOrder(BondOrder::Triple),
            (true, b':') => BondExpression::NotOrder(BondOrder::Aromatic),
            _ => return Err(SmartsError::new(start, "unsupported bond expression")),
        };
        self.pending_bond = Some(expression);
        Ok(())
    }

    fn read_ring(&mut self) -> Result<(), SmartsError> {
        let atom = self
            .current
            .ok_or_else(|| self.error("ring closure has no atom"))?;
        let label = if self.bytes[self.position] == b'%' {
            let start = self.position + 1;
            let end = start + 2;
            let token = self
                .text
                .get(start..end)
                .ok_or_else(|| self.error("ring label after % needs two digits"))?;
            self.position = end;
            token
                .parse::<u16>()
                .map_err(|_| self.error("invalid ring label"))?
        } else {
            let label = u16::from(self.bytes[self.position] - b'0');
            self.position += 1;
            label
        };
        let bond = self.pending_bond.take();
        if let Some((other, opening_bond)) = self.rings.remove(&label) {
            if atom == other {
                return Err(self.error("ring closure forms a self bond"));
            }
            if bond.is_some() && opening_bond.is_some() && bond != opening_bond {
                return Err(self.error("conflicting ring bond expressions"));
            }
            let expression = match bond.or(opening_bond) {
                Some(expression) => expression,
                None => BondExpression::Default,
            };
            self.bonds.push(PatternBond {
                first: other,
                second: atom,
                expression,
            });
        } else {
            self.rings.insert(label, (atom, bond));
        }
        Ok(())
    }

    fn error(&self, message: impl Into<Box<str>>) -> SmartsError {
        SmartsError::new(self.position, message)
    }
}

fn single(test: AtomTest) -> AtomExpression {
    AtomExpression {
        alternatives: vec![vec![SignedAtomTest {
            negated: false,
            test,
        }]],
    }
}

fn parse_bracket(text: &str, offset: usize) -> Result<AtomExpression, SmartsError> {
    if text.is_empty() {
        return Err(SmartsError::new(offset, "empty atom expression"));
    }
    let alternatives = split_top_level(text, b',', offset)?
        .into_iter()
        .map(|(part_offset, part)| parse_conjunction(part, offset + part_offset))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(AtomExpression { alternatives })
}

fn parse_conjunction(
    expression_text: &str,
    offset: usize,
) -> Result<Vec<SignedAtomTest>, SmartsError> {
    let bytes = expression_text.as_bytes();
    let mut position = 0usize;
    let mut tests = Vec::new();
    while position < bytes.len() {
        if matches!(bytes[position], b';' | b'&') {
            position += 1;
            continue;
        }
        let negated = bytes[position] == b'!';
        if negated {
            position += 1;
        }
        let (primitive, consumed) =
            parse_atom_test(&expression_text[position..], offset + position)?;
        tests.push(SignedAtomTest {
            negated,
            test: primitive,
        });
        position += consumed;
    }
    if tests.is_empty() {
        return Err(SmartsError::new(offset, "empty SMARTS alternative"));
    }
    Ok(tests)
}

fn parse_atom_test(text: &str, offset: usize) -> Result<(AtomTest, usize), SmartsError> {
    let bytes = text.as_bytes();
    let first = *bytes
        .first()
        .ok_or_else(|| SmartsError::new(offset, "missing atom test"))?;
    match first {
        b'*' => Ok((AtomTest::Any, 1)),
        b'#' => number_test(&text[1..], offset + 1, |value| {
            AtomTest::Element(Element::from_atomic_number(value))
        }),
        b'a' => Ok((AtomTest::Aromatic, 1)),
        b'A' => Ok((AtomTest::Aliphatic, 1)),
        b'D' => number_test(&text[1..], offset + 1, AtomTest::Degree),
        b'X' => number_test(&text[1..], offset + 1, AtomTest::Connectivity),
        b'H' | b'h' => Ok(optional_number_test(&text[1..], AtomTest::Hydrogens, 1)),
        b'R' => Ok(optional_number_option_test(&text[1..], AtomTest::RingCount)),
        b'r' => Ok(optional_number_option_test(&text[1..], AtomTest::RingSize)),
        b'v' => number_test(&text[1..], offset + 1, AtomTest::Valence),
        b'x' => number_test(&text[1..], offset + 1, AtomTest::RingBonds),
        b'+' | b'-' => parse_charge(text, offset),
        b'@' => Ok(parse_stereo(text)),
        b'$' if bytes.get(1) == Some(&b'(') => parse_recursive(text, offset),
        byte if byte.is_ascii_alphabetic() => parse_element_test(text, offset),
        _ => Err(SmartsError::new(offset, "unsupported atom primitive")),
    }
}

fn parse_element_test(text: &str, offset: usize) -> Result<(AtomTest, usize), SmartsError> {
    let bytes = text.as_bytes();
    let aromatic = bytes[0].is_ascii_lowercase();
    let length =
        if bytes.len() > 1 && bytes[0].is_ascii_uppercase() && bytes[1].is_ascii_lowercase() {
            2
        } else {
            1
        };
    let element = Element::from_symbol(&text[..length])
        .ok_or_else(|| SmartsError::new(offset, "unknown element symbol"))?;
    if aromatic {
        // Element and aromaticity are separate implicit tests in SMARTS. Pack
        // both into a recursive one-atom expression so the public AST stays AND-based.
        let expression = SmartsPattern {
            atoms: vec![AtomExpression {
                alternatives: vec![vec![
                    SignedAtomTest {
                        negated: false,
                        test: AtomTest::Element(element),
                    },
                    SignedAtomTest {
                        negated: false,
                        test: AtomTest::Aromatic,
                    },
                ]],
            }],
            bonds: Vec::new(),
        };
        Ok((AtomTest::Recursive(Box::new(expression)), length))
    } else {
        Ok((AtomTest::Element(element), length))
    }
}

fn number_test(
    text: &str,
    offset: usize,
    constructor: impl FnOnce(u8) -> AtomTest,
) -> Result<(AtomTest, usize), SmartsError> {
    let (value, digits) =
        read_number(text).ok_or_else(|| SmartsError::new(offset, "number required"))?;
    Ok((constructor(value), digits + 1))
}

fn optional_number_test(
    text: &str,
    constructor: impl FnOnce(u8) -> AtomTest,
    default: u8,
) -> (AtomTest, usize) {
    let (value, digits) = match read_number(text) {
        Some(number) => number,
        None => (default, 0),
    };
    (constructor(value), digits + 1)
}

fn optional_number_option_test(
    text: &str,
    constructor: impl FnOnce(Option<u8>) -> AtomTest,
) -> (AtomTest, usize) {
    let number = read_number(text);
    (
        constructor(number.map(|(value, _)| value)),
        1 + number.map_or(0, |(_, digits)| digits),
    )
}

fn read_number(text: &str) -> Option<(u8, usize)> {
    let digits = text.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return None;
    }
    text[..digits].parse().ok().map(|value| (value, digits))
}

fn parse_charge(text: &str, offset: usize) -> Result<(AtomTest, usize), SmartsError> {
    let sign = if text.as_bytes()[0] == b'+' {
        1i8
    } else {
        -1i8
    };
    let repeated = text
        .bytes()
        .take_while(|byte| *byte == text.as_bytes()[0])
        .count();
    let number = read_number(&text[1..]);
    if repeated > 1 && number.is_some() {
        return Err(SmartsError::new(
            offset,
            "mixed repeated and numeric charge",
        ));
    }
    let (magnitude, consumed) = match number {
        Some((value, digits)) => (
            i8::try_from(value).map_err(|_| SmartsError::new(offset, "charge too large"))?,
            1 + digits,
        ),
        None => (
            i8::try_from(repeated).map_err(|_| SmartsError::new(offset, "charge too large"))?,
            repeated,
        ),
    };
    Ok((AtomTest::Charge(sign * magnitude), consumed))
}

fn parse_stereo(text: &str) -> (AtomTest, usize) {
    if text.as_bytes().get(1) == Some(&b'@') {
        (AtomTest::Stereo(StereoConfiguration::S), 2)
    } else {
        (AtomTest::Stereo(StereoConfiguration::R), 1)
    }
}

fn parse_recursive(text: &str, offset: usize) -> Result<(AtomTest, usize), SmartsError> {
    let mut depth = 0usize;
    for (position, byte) in text.bytes().enumerate().skip(1) {
        match byte {
            b'(' => depth += 1,
            b')' => {
                if depth == 0 {
                    return Err(SmartsError::new(
                        offset + position,
                        "unmatched closing parenthesis",
                    ));
                }
                depth -= 1;
                if depth == 0 {
                    let pattern = parse(&text[2..position])?;
                    return Ok((AtomTest::Recursive(Box::new(pattern)), position + 1));
                }
            }
            _ => {}
        }
    }
    Err(SmartsError::new(offset, "unclosed recursive SMARTS"))
}
