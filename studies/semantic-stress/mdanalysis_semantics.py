"""What MDAnalysis does with the PDB-format conditions of the tool-semantics table.

MDAnalysis has no mmCIF reader, so the author/label condition does not apply to it. The
others are read from the same PDB files the other tools read. Run under a Python MDAnalysis
supports:

    uv run --python 3.12 --no-project --with MDAnalysis==2.10.0 \
        python studies/semantic-stress/mdanalysis_semantics.py \
            --refined-set /path/to/refined-set > studies/semantic-stress/mdanalysis_semantics.json
"""

import argparse
import json
import sys
import tempfile
import warnings
from pathlib import Path

import MDAnalysis as mda  # noqa: N813
from build_corpus import verify_files

warnings.simplefilter("ignore")
BENCH = Path("/Users/mcs/Documents/code/biology/molframe/crates/molframe-bench/data")

ALTLOC_PROBE = """\
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C
ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  9.62           C
ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  9.62           O
ATOM      5  CB AMET A   1      25.112  24.880   3.649  0.60 13.77           C
ATOM      6  CB BMET A   1      25.612  24.380   3.049  0.40 13.77           C
ATOM      7  CG AMET A   1      25.353  24.860   5.134  0.60 16.29           C
ATOM      8  CG BMET A   1      26.353  23.860   4.134  0.40 16.29           C
END
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined-set", type=Path, required=True)
    arguments = parser.parse_args()
    corpus = json.loads(Path(__file__).with_name("corpus.json").read_text())
    verify_files(arguments.refined_set, [f for e in corpus["entries"] for f in e["files"].values()])
    verify_files(Path(), corpus["bench"])
    rows = []

    with_icodes = next(e for e in corpus["entries"] if e["facts"]["insertion_code_residues"] >= 5)
    path = arguments.refined_set / with_icodes["files"]["protein"]["path"]
    universe = mda.Universe(str(path))
    residues = universe.residues
    rows.append(
        {
            "condition": f"insertion codes ({with_icodes['id']})",
            "tool": "MDAnalysis",
            "default_behaviour": "a residue is a run of atoms with the same (resid, resname, segid)",
            "information_retained": "icode kept as a per-atom attribute (`icodes`)",
            "result": {
                "residues": len(residues),
                "distinct_resid": len({int(r.resid) for r in residues}),
                "distinct_resid_and_icode": len({(int(r.resid), str(r.icode)) for r in residues}),
            },
        }
    )

    ensemble = mda.Universe(str(BENCH / "2m7c.pdb"))
    rows.append(
        {
            "condition": "32-model NMR ensemble (2M7C)",
            "tool": "MDAnalysis",
            "default_behaviour": "models become trajectory frames; the Universe sits on the first",
            "information_retained": "every model, as frames",
            "result": {"frames": len(ensemble.trajectory), "atoms": len(ensemble.atoms)},
        }
    )

    with tempfile.TemporaryDirectory() as scratch:
        probe = Path(scratch) / "probe.pdb"
        probe.write_text(ALTLOC_PROBE)
        universe = mda.Universe(str(probe))
        rows.append(
            {
                "condition": "alternate conformations (synthetic two-conformer probe built from 1UBQ MET1)",
                "tool": "MDAnalysis",
                "default_behaviour": "every conformer is a separate atom; no conformer is chosen",
                "information_retained": "altLoc per atom (`altLocs`) and occupancy",
                "result": {
                    "atoms": len(universe.atoms),
                    "altlocs": sorted({str(a) for a in universe.atoms.altLocs}),
                },
            }
        )

    unique = Path(arguments.refined_set) / with_icodes["files"]["protein"]["path"]
    hydrogens = mda.Universe(str(unique))
    elements = [str(e).strip().upper() for e in getattr(hydrogens.atoms, "elements", [])]
    names = [str(n).upper() for n in hydrogens.atoms.names]
    rows.append(
        {
            "condition": "explicit hydrogens in the input",
            "tool": "MDAnalysis",
            "default_behaviour": "hydrogens are atoms; selections such as 'not name H*' remove them by name",
            "information_retained": "all atoms and the element column",
            "result": {
                "atoms": len(hydrogens.atoms),
                "hydrogens_by_element_column": int(sum(e == "H" for e in elements)),
                "hydrogens_by_name_prefix": int(sum(n.startswith("H") for n in names)),
            },
        }
    )
    json.dump({"MDAnalysis": mda.__version__, "table": rows}, sys.stdout, indent=1, default=str)
    print()


if __name__ == "__main__":
    main()
