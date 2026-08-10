//! B-factor distribution and explicit TLS consistency validation.

mod distribution;
mod tls;

pub use distribution::{BFactorDistribution, BFactorError, BFactorOutlier, b_factor_distribution};
pub use tls::{TlsBFactorFlag, TlsBFactorReport, TlsGroup, TlsModel, tls_b_factor_consistency};

#[cfg(test)]
#[path = "../bfactor_tests.rs"]
mod tests;
