//! Core-CIF small-molecule structures.

mod model;

pub use model::{
    SmallCifAtom, SmallCifBond, SmallCifDialect, SmallCifError, SmallCifOptions, SmallCifStructure,
    lower_small_cif, lower_small_cif_with_options,
};
