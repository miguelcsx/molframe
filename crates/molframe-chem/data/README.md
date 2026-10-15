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

## Component dictionaries

No component dictionary is compiled in. A dictionary is resolved, in order,
from the structure's own `chem_comp_*` rows, an explicit `--ccd` file, and the
cache entry the profile names with a `components:` line (`<name>.cif`; the file
stem is the recorded version). Fetch the full wwPDB dictionary from
https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz into the cache
with `molframe data fetch` and a manifest, or point `--ccd` at a file. A
component none of these layers defines is reported (W3201), never guessed.

## Embedded default tables

Four files here are byte-identical copies of the repository's `data/` tables
(all MIT), kept beside the crate so the published package carries them. The
opt-in `embedded-tables` feature compiles them in with `include_bytes!` as the
table store's last layer, with a profile naming each; without the feature they
are never read. They are the tables the default flows select when nothing else
does: bond perception on read, bond inference, and surface area and clash
scoring under the default radius policy and metal fallback.

| File | Kind | Bytes | SHA-256 |
| --- | --- | ---: | --- |
| `macromolecular-reach.reach` | `bond-reach` | 664 | `81ed2c2c55bd1b6d98280ef66a376cdb23f4068f5c86260dea24ee51d504190e` |
| `pyykko-covalent.covalent` | `covalent-radii` | 1154 | `733ddda0d7cc7d8698d688c39fdf05a75ead78d78e79221810ce3a260f4a8710` |
| `bondi.radii` | `vdw-radii` | 413 | `6ac4e602e667b62ce7597f75b1367c38fbcecf0d53c4680be9edefd00c22d844` |
| `macromolecular-metals.metals` | `metals` | 636 | `c1b08c545099c6368fc0baeb2127bf1463cbaf685e80efd1a143bdfd99ec865f` |

`src/tables/table_bundled_tests.rs` keeps every copy equal to its source and
every compiled-in entry equal to the file it names. Edit the root `data/` file,
then copy it here.
