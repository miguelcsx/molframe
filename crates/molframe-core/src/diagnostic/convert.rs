//! Turning a kernel's own error into a registered diagnostic.
//!
//! Every kernel keeps a small error enum of its own, and a caller that wants one
//! vocabulary for all of them — the Python bindings above all — needs each to
//! name a stable code. The mapping lives beside the enum so a new variant must
//! choose its code where it is declared.

/// Implements `From<&E>` and `From<E>` for [`Diagnostic`](crate::Diagnostic).
///
/// The code is chosen by the expression given and the message is the error's
/// own `Display` text, so nothing the error already says is lost. An optional
/// third argument adds structured context to the diagnostic.
///
/// ```
/// use molframe_core::{Code, Diagnostic, diagnostic_from};
///
/// #[derive(Debug)]
/// struct TooBig(usize);
/// impl std::fmt::Display for TooBig {
///     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
///         write!(f, "{} is too big", self.0)
///     }
/// }
/// diagnostic_from!(TooBig, |_error| Code::E7001, |diagnostic, error| {
///     diagnostic.with_context("size", error.0.to_string())
/// });
///
/// let diagnostic = Diagnostic::from(TooBig(7));
/// assert_eq!(diagnostic.code(), Code::E7001);
/// ```
#[macro_export]
macro_rules! diagnostic_from {
    ($error:ty, |$value:ident| $code:expr) => {
        $crate::diagnostic_from!($error, |$value| $code, |diagnostic, _error| diagnostic);
    };
    ($error:ty, |$value:ident| $code:expr, |$diagnostic:ident, $inner:ident| $context:expr) => {
        impl ::core::convert::From<&$error> for $crate::Diagnostic {
            fn from($value: &$error) -> Self {
                let code: $crate::Code = $code;
                let $diagnostic = $crate::Diagnostic::new(code)
                    .with_message(::std::string::ToString::to_string($value));
                let $inner = $value;
                $context
            }
        }

        impl ::core::convert::From<$error> for $crate::Diagnostic {
            fn from(value: $error) -> Self {
                Self::from(&value)
            }
        }
    };
}

impl From<&std::convert::Infallible> for crate::Diagnostic {
    /// An error that cannot be constructed has no code to name.
    fn from(never: &std::convert::Infallible) -> Self {
        match *never {}
    }
}

impl From<&crate::Diagnostic> for crate::Diagnostic {
    /// A finding is already a diagnostic; kernels that fail with one convert by copy.
    fn from(diagnostic: &crate::Diagnostic) -> Self {
        diagnostic.clone()
    }
}
