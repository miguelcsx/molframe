# Data

Small permissive tables that molframe reads at run time. None is compiled
into a binary; each is a file in the shared text format (a `# SPDX` header,
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

`build_reference.py` regenerates the Top8000 tables from the upstream files.
