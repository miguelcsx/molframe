# Element reference data

The generated tables pin mendeleev 1.1.0 as the distribution source. The
scientific series retain their own identities in the Rust API:

- single-bond covalent radii: Pyykkö and Atsumi (2009);
- Pauling electronegativity values: CRC compilation;
- van der Waals radii: Bondi (1964) and Alvarez (2013);
- ionic and crystal radii: Shannon (1976), including charge, coordination and
  spin state.

Missing source values remain `None`; one convention is never substituted for
another. Regenerate from the pinned `elements.db` with
`scripts/generate-elements.py`. Mendeleev's MIT license is retained in
`LICENSE.mendeleev`.
