# Data

Small permissive tables that molframe reads at run time. Only one is ever
compiled into a binary: `macromolecular-reach.reach`, copied to
`crates/molframe-chem/data/` and compiled in by the opt-in `embedded-tables`
feature (a test keeps the copies identical). Each is a file in the shared text format (a `# SPDX` header,
`format: 1`, `name`, `version`, `source`, `licence`, then rows) or the binary
form of it (`MFR2`/`MFT2` magic, a header, a payload).

| File | Kind | Licence | Origin |
| --- | --- | --- | --- |
| `top8000-rama.rama` | ramachandran | CC-BY-4.0 | Top8000, Richardson Laboratory, Duke University |
| `top8000-rotamer.rotamer` | rotamer | CC-BY-4.0 | Top8000, Richardson Laboratory, Duke University |
| `alphafold-peptide-link.geometry` | geometry | Apache-2.0 | AlphaFold `residue_constants.py`, peptide-link terms only |
| `chi-atoms.chi` | chi | MIT | side-chain chi atom definitions |
| `nucleic-sanity.nucleic` | nucleic | MIT | sugar-pucker and backbone windows |
| `heavy-atom-roles.roles` | heavy-roles | MIT | donor and acceptor roles by residue and atom name |
| `tien-2013-max-area.maxarea` | max-area | CC-BY-4.0 | Tien et al. 2013 maximum accessible areas |
| `bondi.radii` | vdw-radii | MIT | Bondi (1964), via mendeleev 1.1.0 |
| `alvarez.radii` | vdw-radii | MIT | Alvarez (2013) van der Waals radii |
| `charmm.radii`, `amber-united.radii` | vdw-radii | MIT | CHARMM36m and AMBER united-atom radii by element |
| `pyykko-covalent.covalent` | covalent-radii | MIT | Pyykko and Atsumi (2009), via mendeleev 1.1.0 |
| `shannon-ionic.ionic` | ionic-radii | MIT | Shannon (1976), via mendeleev 1.1.0 |
| `macromolecular-metals.metals` | metals | MIT | typical metal charge and coordination in macromolecules; radii from Alvarez |
| `macromolecular-reach.reach` | bond-reach | MIT | per-element reach for inferred bonds |

## Attribution

The Ramachandran and rotamer tables derive from the Top8000 reference data of
the Richardson Laboratory, Duke University, licensed CC-BY-4.0 (see
`LICENSES/CC-BY-4.0.txt`). molframe reformatted them into the table format;
the licence, source and version travel in each table's identity and appear in
every result that applied it.

## Use

Copy the files into the cache directory (`MOLFRAME_DATA_DIR`, else
`$XDG_CACHE_HOME/molframe`, else `~/.cache/molframe`) and list the entry per
kind in `profile.txt` (`kind: name`). Tables that cannot be redistributed, such
as the full bond and angle geometry, come from a manifest:
`molframe data fetch <entry> --manifest <file>`.

A table is selected by name from the profile (`kind: name`, kinds `vdw-radii`,
`covalent-radii`, `ionic-radii`, `metals`, `bond-reach`) or by passing a file
path, or registered in memory with `TableStore::insert_profile` and
`insert_table`. There is no built-in default except the bonding reach under
`embedded-tables`, which a registration or a profile on disk overrides; a
missing table is an `E7102` error.
`crates/molframe-bench/data/profile.txt` is the profile the tests and benches
use: it names the tables in this directory, and `.cargo/config.toml` points
`MOLFRAME_DATA_DIR` here and `MOLFRAME_PROFILE` at it.

`build_reference.py` regenerates the Top8000 tables from the upstream files.
