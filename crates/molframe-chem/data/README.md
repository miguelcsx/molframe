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

## Standard amino-acid dictionary fixture

`CCD-amino-acids.cif` concatenates unmodified wwPDB CCD records for the twenty
standard amino acids, retrieved from the RCSB ligand endpoint
(https://files.rcsb.org/ligands/download/{ID}.cif) on 2026-10-03. SHA-256:
`7cc6f3143945986ad1f0e8913f8e5980f966e14ec44a6bfed166986b9bc6db0f`.
The data are CC0-1.0 (https://www.wwpdb.org/about/usage-policies). It is a
protein-only dictionary: structures with heme, waters or other hetero residues
need those components supplied too, and the comparison and mapping code refuses
a residue the dictionary does not define rather than guessing.
