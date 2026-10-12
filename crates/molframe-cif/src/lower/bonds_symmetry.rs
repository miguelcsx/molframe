//! Whether a `struct_conn` partner sits on the asymmetric unit or on an image.
//!
//! A connection may name its partner through a symmetry operator and a
//! translation of the unit cell. The bond graph holds atoms of one asymmetric
//! unit and has no way to name an image, so a connection whose partner is
//! genuinely on a neighbouring image must not be lowered to a bond between the
//! local atoms that happen to share those identifiers — that would assert a
//! contact the deposition does not describe.
//!
//! The test here is deliberately narrow: it accepts a partner the block shows to
//! be the local atom and refuses everything else. It does not do crystallography;
//! composing operators and placing images belongs to a reader that carries cells
//! and operators.

use crate::document::DataBlock;

/// The symmetry operators the block defines, in file order.
///
/// Both the current and the legacy category are read; the current one wins when a
/// block carries both.
pub(super) struct OperatorTable {
    operations: Vec<Box<str>>,
}

impl OperatorTable {
    /// Collects the operators the block defines, if it defines any.
    pub(super) fn of(block: &DataBlock) -> Self {
        let category = block
            .category("space_group_symop")
            .map(|category| (category, "operation_xyz"))
            .or_else(|| {
                block
                    .category("symmetry_equiv")
                    .map(|category| (category, "pos_as_xyz"))
            });
        let mut operations = Vec::new();
        if let Some((category, item)) = category {
            for row in 0..category.row_count() {
                if let Some(text) = category.text(item, row) {
                    operations.push(Box::from(text));
                }
            }
        }
        Self { operations }
    }

    /// The operator numbered `number`, counting from one.
    fn operation(&self, number: u32) -> Option<&str> {
        let position = usize::try_from(number).ok()?.checked_sub(1)?;
        self.operations.get(position).map(Box::as_ref)
    }

    fn is_empty(&self) -> bool {
        self.operations.is_empty()
    }
}

/// A partner that cannot be expressed as an atom of this asymmetric unit.
pub(super) struct Unrepresentable<'a> {
    /// Which partner carried the value.
    pub(super) partner: &'static str,
    /// The value as the file wrote it.
    pub(super) value: &'a str,
}

/// The partner and value that make a connection unrepresentable, if any.
///
/// An absent, unknown or inapplicable value says nothing about the partner and is
/// read as local data the caller already had. Anything explicit has to be shown to
/// be the local atom.
pub(super) fn unrepresentable<'a>(
    first: Option<&'a str>,
    second: Option<&'a str>,
    operators: &OperatorTable,
) -> Option<Unrepresentable<'a>> {
    for (partner, value) in [("ptnr1_symmetry", first), ("ptnr2_symmetry", second)] {
        let Some(value) = value.filter(|text| is_stated(text)) else {
            continue;
        };
        if !is_local(value, operators) {
            return Some(Unrepresentable { partner, value });
        }
    }
    None
}

/// Whether a field states a value rather than marking it unknown.
///
/// The dictionary's `?` and `.` are not symmetry codes; they are the absence of
/// one, which the deposited archive uses for a partner already inside the
/// asymmetric unit.
fn is_stated(value: &str) -> bool {
    let trimmed = value.trim();
    !trimmed.is_empty() && trimmed != "?" && trimmed != "."
}

/// Whether an explicit symmetry value names the local asymmetric unit.
///
/// The code is `n_klm` or a bare operator number `n`. The translation digits are
/// stored offset by five, so `n_555` is no translation at all. An operator other
/// than the first is never the identity on its own, and with an operator table in
/// hand the named row has to say so in its own right — a table is free to define
/// operator 1 as something other than `x,y,z`.
fn is_local(value: &str, operators: &OperatorTable) -> bool {
    let value = value.trim();
    let Some((operator, translation)) = parse(value) else {
        return false;
    };
    if translation != [0, 0, 0] {
        return false;
    }
    if operators.is_empty() {
        // No table to consult: operator one is the identity by the archive's own
        // convention, and nothing else can be assumed.
        return operator == 1;
    }
    operators
        .operation(operator)
        .is_some_and(is_identity_operation)
}

/// Splits `n_klm` into the operator number and the cell translation it applies.
///
/// The translation is returned in unit-cell repeats, so `_555` is `[0, 0, 0]` and
/// `_565` is `[0, 1, 0]`.
fn parse(value: &str) -> Option<(u32, [i32; 3])> {
    match value.split_once('_') {
        Some((operator, digits)) => {
            if digits.len() != 3 {
                return None;
            }
            let mut translation = [0i32; 3];
            for (slot, character) in translation.iter_mut().zip(digits.bytes()) {
                // The archive offset digits run from zero to nine, standing for -5
                // to +4. Anything else is not a code this test can read.
                if !character.is_ascii_digit() {
                    return None;
                }
                let digit = u32::from(character) - u32::from(b'0');
                *slot = i32::try_from(digit).ok()?.checked_sub(5)?;
            }
            Some((operator.parse().ok()?, translation))
        }
        None => Some((value.parse().ok()?, [0, 0, 0])),
    }
}

/// Whether an operator expression places a point where it already was.
///
/// Only the literal identity coordinates count. Evaluating an arbitrary
/// expression would mean implementing crystallography here, and a wrong "yes"
/// invents a bond, so anything else is refused.
fn is_identity_operation(operation: &str) -> bool {
    // The identity expression is ten bytes at most, so a longer one cannot match
    // and the fixed buffer keeps this off the heap on a per-connection path.
    let mut buffer = [0u8; 10];
    let mut length = 0usize;
    for byte in operation.bytes() {
        if byte.is_ascii_whitespace() {
            continue;
        }
        let Some(slot) = buffer.get_mut(length) else {
            return false;
        };
        *slot = byte.to_ascii_lowercase();
        length += 1;
    }
    let Some(text) = buffer.get(..length) else {
        return false;
    };
    let text = match text.strip_prefix(b"+".as_slice()) {
        Some(rest) => rest,
        None => text,
    };
    matches!(text, b"x,y,z" | b"x,y,z,1" | b"x,y,z,1.0")
}

#[cfg(test)]
#[path = "bonds_symmetry_tests.rs"]
mod tests;
