//! Substitution matrices as reference tables, and the dense lookup built from one.

mod load;
mod score_matrix;
mod table;

pub use load::load_matrix;
pub use score_matrix::ScoreMatrix;
pub use table::SubstitutionMatrix;
