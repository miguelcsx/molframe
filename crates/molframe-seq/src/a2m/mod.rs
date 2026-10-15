//! A2M alignment format.

mod codec;

pub use codec::{A2mError, a2m_match_columns, parse_a2m, write_a2m};
