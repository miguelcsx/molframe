//! Process exit codes.
//!
//! The codes are stable and documented, because a script must be able to tell
//! "the answer is X" from "there is no defensible answer" without reading prose.
//! That distinction is the whole point of separating a refused conversion from a
//! generic failure, and an indeterminate analysis from either.

use pdbiox::{Code, Diagnostic};

/// What the process exits with.
///
/// The whole set is defined here rather than only the codes something currently
/// produces, because the numbers are a published contract: a script branches on
/// them, and a code that appears later must not shift the ones already in use.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(i32)]
pub enum Exit {
    /// Everything worked.
    Success = 0,
    /// Something failed that has no more specific code.
    Failure = 1,
    /// The arguments were wrong.
    Usage = 2,
    /// The input could not be found or read.
    Input = 3,
    /// The input could not be parsed.
    Parse = 4,
    /// The input violated the schema it claims to follow.
    Schema = 5,
    /// The structure is internally inconsistent.
    Consistency = 6,
    /// The conversion was refused because data would have been lost.
    Refused = 7,
    /// The policy was contradictory or named something that does not exist.
    Policy = 8,
    /// The analysis declined to produce an answer.
    Indeterminate = 9,
    /// A resource limit was exceeded.
    Resource = 10,
}

impl Exit {
    /// The code a set of findings implies.
    ///
    /// Conversion refusals outrank everything else, because a refused write is
    /// the specific thing a caller most needs to distinguish, and a run that
    /// refused is not a run that merely failed.
    #[must_use]
    pub fn of(findings: &[Diagnostic]) -> Self {
        let mut worst = Self::Success;
        for finding in findings {
            let candidate = Self::of_code(finding.code());
            if candidate == Self::Refused {
                return Self::Refused;
            }
            if candidate as i32 > worst as i32 {
                worst = candidate;
            }
        }
        worst
    }

    fn of_code(code: Code) -> Self {
        use pdbiox::Class;
        match code.class() {
            Class::Syntax => match code {
                Code::E1901 => Self::Resource,
                Code::E1001 => Self::Input,
                _ => Self::Parse,
            },
            Class::Schema => Self::Schema,
            Class::Consistency => Self::Consistency,
            Class::Conversion => Self::Refused,
            Class::Policy => Self::Policy,
            Class::Resource => Self::Resource,
            // A class this build does not know about is still a failure, and
            // saying so beats refusing to compile against a newer library.
            Class::Geometry | Class::Internal | _ => Self::Failure,
        }
    }

    /// The value to hand to the operating system.
    #[must_use]
    pub const fn code(self) -> i32 {
        self as i32
    }
}

#[cfg(test)]
#[path = "exit_tests.rs"]
mod tests;
