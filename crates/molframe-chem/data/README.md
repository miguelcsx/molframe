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

## Carbohydrate reference metadata and fixtures

`src/carbohydrates/snfg.rs` adapts all 79 Mol* `Monosaccharides` metadata records
and its `CommonSaccharideNames` CCD aliases, pinned to Mol* commit
`e58afe3353b9d6346b78ec5ad9f36fe0e4cd8789`, file
`src/mol-model/structure/structure/carbohydrates/constants.ts` (SHA-256
`90f0f7de929026248920f50a08c1e68827a23c59bce4a17f56467237beca941e`).
Copyright 2018–2026 Mol* contributors; MIT permission notice: `LICENSE.molstar`.
Shapes/colours follow https://www.ncbi.nlm.nih.gov/glycans/snfg.html.
CHARMM/GLYCAM aliases are excluded: their namespaces collide with CCD codes
(for example 4GL and UEA). Unknown sugars use Mol* neutral flat hexagons.

The complete deposited 1HZH CIF is stored losslessly as `1HZH.cif.gz` from
https://files.rcsb.org/download/1HZH.cif, DOI 10.2210/pdb1HZH/pdb. Decompressed
SHA-256: `5793921cc7dcc2317e1892edcbf32b2a8e11a7bef81d2eb25d94c6165f0b01ab`;
gzip SHA-256: `933296767ba9eb5b2ff2b930c9501af301ec6ec116f7d4c704a1f7a9355dc196`.
`CCD-saccharides.cif` concatenates unmodified NAG/MAN/BMA/FUC/GAL/ASN CCD files
from https://files.rcsb.org/ligands/download/{ID}.cif, acquired 2026-10-02;
SHA-256: `ff87cfab53b91d47ee64e3782df4f8d06b91c1278ac5e16c01ec99c9214e5fc4`.
ASN supplies native acceptor topology for conservative attachment inference;
no glycan–protein attachment is declared by 1HZH `struct_conn`.
These wwPDB data are CC0-1.0 (https://www.wwpdb.org/about/usage-policies).

## Peptide charge fixture

`CCD-amino-acids.cif` concatenates unmodified wwPDB CCD records for the twenty
standard amino acids, retrieved from the same RCSB ligand endpoint on
2026-10-03. SHA-256:
`7cc6f3143945986ad1f0e8913f8e5980f966e14ec44a6bfed166986b9bc6db0f`.
The data are CC0-1.0. Charge validation uses the committed 4HHB protein atoms
with an explicit dictionary; heme, waters and unparameterised metal chemistry
are excluded from this protein-only fixture, not silently assigned zero.
