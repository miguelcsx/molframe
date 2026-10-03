//! Topological saccharide recognition and conservative glycosidic connectivity.
//!
//! Rings contain one oxygen and four or five carbons. Recognition never infers
//! ring bonds from proximity. The optional linkage fallback reuses the chemistry
//! cell grid and considers only vacant anomeric carbons and chemically eligible
//! acceptors within 2 Å. CSR construction is O(N+B); bounded-depth ring traversal
//! is O(S d^5), with chemically bounded degree d, and the grid is O(N k) at bounded
//! local density. Storage is O(N+B+S+L).

mod geometry;
mod links;
mod perceive;
mod rings;
mod snfg;
mod types;

pub use perceive::carbohydrates;
pub use snfg::snfg_symbol;
pub use types::{
    CarbohydrateLink, CarbohydrateOptions, CarbohydrateReport, Monosaccharide, RingGeometry,
    SnfgShape, SnfgSymbol,
};
