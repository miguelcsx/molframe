# Element reference data

The tables are pinned to mendeleev 1.1.0 as the distribution source, and are
committed as `crates/molframe-chem/src/element_data.rs` — a hand-maintained data
module, not a build step. The scientific series retain their own identities in
the Rust API:

- single-bond covalent radii: Pyykkö and Atsumi (2009);
- Pauling electronegativity values: CRC compilation;
- van der Waals radii: Bondi (1964) and Alvarez (2013);
- ionic and crystal radii: Shannon (1976), including charge, coordination and
  spin state.

Missing source values remain `None`; one convention is never substituted for
another. Mendeleev's MIT license is retained in `LICENSE.mendeleev`.
