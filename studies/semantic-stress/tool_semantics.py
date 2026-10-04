"""What each tool does with the same structural condition, observed rather than recalled.

For each condition in the corpus this reads the same file with MolFrame,
Biopython, Biotite and Gemmi (and MDAnalysis where it reads the format), asks
each the same question, and records what the tool exposes by default, what it
keeps reachable, and the number it gives. Nothing here judges a tool: a default
is a choice, and the point is that the choices differ.

Run from the repository root:

    PYTHONPATH=python uv run --no-sync \
        --with biopython==1.88 --with biotite==1.7.1 --with gemmi==0.7.5 \
        python studies/semantic-stress/tool_semantics.py \
            --refined-set /path/to/refined-set > studies/semantic-stress/tool_semantics.json
"""

import argparse
import gzip
import io
import json
import sys
import tempfile
import warnings
from collections.abc import Callable
from importlib import metadata
from pathlib import Path

import numpy as np
from build_corpus import verify_files

import molframe

warnings.simplefilter("ignore")
BENCH = Path("crates/molframe-bench/data")


def versions() -> dict[str, str]:
    found = {"molframe": molframe.__version__ if hasattr(molframe, "__version__") else "workspace"}
    for package in ("biopython", "biotite", "gemmi", "MDAnalysis"):
        try:
            found[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            found[package] = "not installed"
    return found


def row(case: str, tool: str, default: str, retained: str, result: object) -> dict[str, object]:
    return {
        "condition": case,
        "tool": tool,
        "default_behaviour": default,
        "information_retained": retained,
        "result": result,
    }


def guarded(case: str, tool: str, run: Callable[[], dict[str, object]]) -> dict[str, object]:
    try:
        return row(case, tool, **run())
    except Exception as error:  # noqa: BLE001 - a tool failing on an input is itself the observation
        return row(case, tool, "raised", "n/a", f"{type(error).__name__}: {error}")


# ---------------------------------------------------------------------------------------------
# author versus label identifiers (1AON mmCIF)
# ---------------------------------------------------------------------------------------------


def divergent_residue() -> dict[str, object]:
    """The first residue of 1AON whose author and label sequence numbers differ."""
    structure = molframe.read(BENCH / "1aon.cif.gz")
    for chain_index in range(len(structure.models[0].chains)):
        chain = structure.models[0].chains[chain_index]
        for position in range(len(chain.residues)):
            residue = chain.residues[position]
            if residue.number is not None and residue.auth_number not in (None, residue.number):
                return {
                    "label_chain": chain.label,
                    "auth_chain": chain.auth_label,
                    "label_seq": residue.number,
                    "auth_seq": residue.auth_number,
                    "name": residue.name,
                }
    msg = "no divergent residue found"
    raise RuntimeError(msg)


def case_identifiers() -> list[dict[str, object]]:
    case = "auth_seq_id differs from label_seq_id (1AON)"
    target = divergent_residue()
    path = BENCH / "1aon.cif.gz"
    text = gzip.decompress(path.read_bytes()).decode()
    rows: list[dict[str, object]] = []

    def molframe_run() -> dict[str, object]:
        structure = molframe.read(path)
        by_auth = structure.select(f"chain {target['auth_chain']} and resid {target['auth_seq']}")
        by_label = structure.select(
            f"label_chain {target['label_chain']} and label_resid {target['label_seq']}"
        )
        return {
            "default": "selectors read author ids (Namespace::Auth); label ids need an explicit selector",
            "retained": "both, as separate columns, per residue",
            "result": {
                "atoms_by_author_ids": len(by_auth.indices),
                "atoms_by_label_ids": len(by_label.indices),
            },
        }

    def biopython_run() -> dict[str, object]:
        from Bio.PDB.MMCIFParser import MMCIFParser

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("x", io.StringIO(text))
        chain = structure[0][str(target["auth_chain"])]
        found_auth = [r for r in chain if r.id[1] == target["auth_seq"]]
        found_label = [r for r in chain if r.id[1] == target["label_seq"]]
        return {
            "default": "residue.id[1] is auth_seq_id; chain.id is auth_asym_id",
            "retained": "label ids only in the raw MMCIF2Dict, not on Structure objects",
            "result": {
                "atoms_by_author_ids": sum(len(r) for r in found_auth),
                "atoms_by_label_ids": sum(len(r) for r in found_label),
                "same_number_names_a_different_residue": bool(found_label)
                and found_label != found_auth,
            },
        }

    def biotite_run() -> dict[str, object]:
        import biotite.structure.io.pdbx as pdbx

        cif = pdbx.CIFFile.read(io.StringIO(text))
        author = pdbx.get_structure(cif, model=1, use_author_fields=True)
        label = pdbx.get_structure(cif, model=1, use_author_fields=False)

        def count(atoms, chain, seq):
            return int(((atoms.chain_id == chain) & (atoms.res_id == seq)).sum())

        return {
            "default": "use_author_fields=True: chain_id and res_id are auth_asym_id and auth_seq_id",
            "retained": "either, by re-reading with use_author_fields=False; not both in one array",
            "result": {
                "atoms_by_author_ids": count(author, target["auth_chain"], target["auth_seq"]),
                "atoms_by_label_ids": count(label, target["label_chain"], target["label_seq"]),
            },
        }

    def gemmi_run() -> dict[str, object]:
        import gemmi

        structure = gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block())
        chain = structure[0][str(target["auth_chain"])]
        found_auth = [r for r in chain if r.seqid.num == target["auth_seq"]]
        found_label = [r for r in chain if r.label_seq == target["label_seq"]]
        return {
            "default": "seqid is auth_seq_id and chain.name is auth_asym_id",
            "retained": "label_seq on every residue, label chain as the residue's subchain",
            "result": {
                "atoms_by_author_ids": sum(len(r) for r in found_auth),
                "atoms_by_label_ids": sum(len(r) for r in found_label),
            },
        }

    rows.append(guarded(case, "molframe", molframe_run))
    rows.append(guarded(case, "biopython", biopython_run))
    rows.append(guarded(case, "biotite", biotite_run))
    rows.append(guarded(case, "gemmi", gemmi_run))
    rows.append(
        row(
            case,
            "MDAnalysis",
            "reads PDB, GRO and other formats; has no mmCIF reader",
            "n/a",
            "not applicable to an mmCIF-only distinction",
        )
    )
    for each in rows:
        each["probe"] = target
    return rows


# ---------------------------------------------------------------------------------------------
# insertion codes (refined set)
# ---------------------------------------------------------------------------------------------


def entry_with_insertions(refined: Path) -> str:
    corpus = json.loads(Path(__file__).with_name("corpus.json").read_text())
    for entry in corpus["entries"]:
        if entry["facts"]["insertion_code_residues"] >= 5:
            return entry["id"]
    msg = "no corpus entry with insertion codes"
    raise RuntimeError(msg)


def case_insertion_codes(refined: Path) -> list[dict[str, object]]:
    code = entry_with_insertions(refined)
    case = f"insertion codes ({code}): residues numbered 52, 52A, 52B"
    path = refined / code / f"{code}_protein.pdb"
    lines = [l for l in path.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    keys = {(l[21], l[22:26].strip(), l[26]) for l in lines}
    true_residues = len(keys)
    by_number = len({(chain, number) for chain, number, _ in keys})
    probe = {"entry": code, "residues": true_residues, "distinct_chain_and_number": by_number}
    rows = []

    def molframe_run() -> dict[str, object]:
        structure = molframe.read(path)
        return {
            "default": "residue identity is the topology row; the insertion code is a column",
            "retained": "insertion code per residue",
            "result": {"residues": structure.residue_count},
        }

    def biopython_run() -> dict[str, object]:
        from Bio.PDB import PDBParser

        structure = PDBParser(QUIET=True).get_structure("x", str(path))
        count = sum(1 for _ in structure[0].get_residues())
        return {
            "default": "residue id is the triple (hetero flag, number, insertion code)",
            "retained": "insertion code in residue.id[2]",
            "result": {"residues": count},
        }

    def biotite_run() -> dict[str, object]:
        import biotite.structure as bs
        import biotite.structure.io.pdb as pdb

        atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
        by_triple = len(bs.get_residue_starts(atoms))
        return {
            "default": "res_id is the number; ins_code is a separate annotation",
            "retained": "ins_code annotation; residue starts are computed from changes of (chain, res_id, ins_code, res_name)",
            "result": {
                "residues": int(by_triple),
                "distinct_chain_and_res_id": len(
                    set(zip(atoms.chain_id, atoms.res_id, strict=True))
                ),
            },
        }

    def gemmi_run() -> dict[str, object]:
        import gemmi

        structure = gemmi.read_structure(str(path))
        count = sum(len(chain) for chain in structure[0])
        return {
            "default": "ResidueId is (number, insertion code)",
            "retained": "icode on every residue",
            "result": {"residues": count},
        }

    rows.append(guarded(case, "molframe", molframe_run))
    rows.append(guarded(case, "biopython", biopython_run))
    rows.append(guarded(case, "biotite", biotite_run))
    rows.append(guarded(case, "gemmi", gemmi_run))
    for each in rows:
        each["probe"] = probe
    return rows


# ---------------------------------------------------------------------------------------------
# explicit hydrogens (refined set): what a surface area depends on
# ---------------------------------------------------------------------------------------------


def case_hydrogens(refined: Path) -> list[dict[str, object]]:
    corpus = json.loads(Path(__file__).with_name("corpus.json").read_text())
    code = sorted(corpus["entries"], key=lambda e: e["facts"]["atoms"])[
        len(corpus["entries"]) // 4
    ]["id"]
    case = f"explicit hydrogens in the input ({code}): heavy-atom SASA with and without them"
    path = refined / code / f"{code}_protein.pdb"
    rows = []

    def molframe_run() -> dict[str, object]:
        structure = molframe.read(path)
        elements = np.asarray(structure.elements)
        radii = np.asarray(molframe.chemistry.vdw_radii(structure), dtype=np.float32)
        xyz = np.asarray(structure.coordinates, dtype=np.float32)
        heavy = elements != 1
        with_h = np.asarray(molframe.surface.sasa(xyz, radii, probe=1.4, points=960))[heavy]
        without = np.asarray(molframe.surface.sasa(xyz[heavy], radii[heavy], probe=1.4, points=960))
        return {
            "default": "every atom of the file occludes; hydrogens are atoms like the rest",
            "retained": "all atoms; removal is an explicit selection",
            "result": {
                "hydrogens": int((elements == 1).sum()),
                "heavy_sasa_with_hydrogens": float(with_h.sum()),
                "heavy_sasa_without_hydrogens": float(without.sum()),
            },
        }

    def biotite_run() -> dict[str, object]:
        import biotite.structure as bs
        import biotite.structure.io.pdb as pdb

        atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
        heavy = atoms[atoms.element != "H"]
        with_h = bs.sasa(atoms, probe_radius=1.4, point_number=960)
        without = bs.sasa(heavy, probe_radius=1.4, point_number=960)
        keep = np.flatnonzero(atoms.element != "H")
        return {
            "default": "every atom occludes, with radii from its own table",
            "retained": "all atoms; filtering is the caller's mask",
            "result": {
                "hydrogens": int((atoms.element == "H").sum()),
                "heavy_sasa_with_hydrogens": float(np.nansum(with_h[keep])),
                "heavy_sasa_without_hydrogens": float(np.nansum(without)),
            },
        }

    rows.append(guarded(case, "molframe", molframe_run))
    rows.append(guarded(case, "biotite", biotite_run))
    for name, why in (
        (
            "biopython",
            "ShrakeRupley uses a radius table keyed by element; hydrogens occlude unless removed",
        ),
        ("gemmi", "has no SASA; neighbour searches include hydrogens unless filtered"),
    ):
        rows.append(row(case, name, why, "all atoms", "not run: no comparable SASA default"))
    return rows


# ---------------------------------------------------------------------------------------------
# ensembles (2M7C)
# ---------------------------------------------------------------------------------------------


def case_ensemble() -> list[dict[str, object]]:
    path = BENCH / "2m7c.pdb"
    case = "32-model NMR ensemble (2M7C)"
    rows = []

    def molframe_run() -> dict[str, object]:
        structure = molframe.read(path)
        return {
            "default": "analyses run on the first model unless the policy names another or all",
            "retained": "every model",
            "result": {"models": len(structure.models), "atoms_per_model": structure.atom_count},
        }

    def biopython_run() -> dict[str, object]:
        from Bio.PDB import PDBParser

        structure = PDBParser(QUIET=True).get_structure("x", str(path))
        return {
            "default": "iterating a Structure yields every model; most examples index model 0",
            "retained": "every model",
            "result": {
                "models": len(structure),
                "atoms_in_model_0": len(list(structure[0].get_atoms())),
            },
        }

    def biotite_run() -> dict[str, object]:
        import biotite.structure.io.pdb as pdb

        stack = pdb.PDBFile.read(str(path)).get_structure()
        return {
            "default": "model=None returns an AtomArrayStack of every model",
            "retained": "every model, as a stack with shared annotations",
            "result": {
                "models": int(stack.stack_depth()),
                "atoms_per_model": int(stack.array_length()),
            },
        }

    def gemmi_run() -> dict[str, object]:
        import gemmi

        structure = gemmi.read_structure(str(path))
        return {
            "default": "Structure holds every model",
            "retained": "every model",
            "result": {
                "models": len(structure),
                "atoms_in_model_0": structure[0].count_atom_sites(),
            },
        }

    rows.append(guarded(case, "molframe", molframe_run))
    rows.append(guarded(case, "biopython", biopython_run))
    rows.append(guarded(case, "biotite", biotite_run))
    rows.append(guarded(case, "gemmi", gemmi_run))
    return rows


# ---------------------------------------------------------------------------------------------
# alternate locations: a synthetic probe, because no corpus entry here has any
# ---------------------------------------------------------------------------------------------


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


def case_altloc() -> list[dict[str, object]]:
    case = "alternate conformations (synthetic two-conformer probe built from 1UBQ MET1)"
    directory = Path(tempfile.mkdtemp())
    path = directory / "probe.pdb"
    path.write_text(ALTLOC_PROBE)
    rows = []

    def molframe_run() -> dict[str, object]:
        from molframe import analysis

        structure = molframe.read(path)
        default = analysis.contacts(structure, 5.0)
        every = analysis.contacts(structure, 5.0, policy=molframe.AnalysisPolicy(altloc="keep-all"))
        return {
            "default": "ConformerConsistent: one self-consistent label per region, applied by every governed analysis and recorded in its coverage",
            "retained": "every conformer as atoms with an altloc label and occupancy; selection alone does not resolve them",
            "result": {
                "atoms_in_file": structure.atom_count,
                "atoms_analysed_by_default": default.coverage.used,
                "status_by_default": default.status,
                "atoms_analysed_keep_all": every.coverage.used,
                "status_keep_all": every.status,
                "atoms_selected_by_a_query_under_the_default": len(structure.select("all").indices),
            },
        }

    def biopython_run() -> dict[str, object]:
        from Bio.PDB import PDBParser

        structure = PDBParser(QUIET=True).get_structure("x", str(path))
        residue = next(iter(structure[0].get_residues()))
        beta = residue["CB"]
        return {
            "default": "disordered atoms are DisorderedAtom; the highest-occupancy conformer is selected",
            "retained": "all conformers inside DisorderedAtom",
            "result": {
                "atoms_iterated": len(list(structure[0].get_atoms())),
                "cb_altloc_selected": beta.get_altloc(),
                "cb_x": round(float(beta.coord[0]), 3),
            },
        }

    def biotite_run() -> dict[str, object]:
        import biotite.structure.io.pdb as pdb

        file = pdb.PDBFile.read(str(path))
        first = file.get_structure(model=1)
        all_ = file.get_structure(model=1, altloc="all")
        return {
            "default": "altloc='first': the first label of each atom",
            "retained": "with altloc='all' every conformer, labelled in an annotation",
            "result": {
                "atoms_by_default": int(first.array_length()),
                "atoms_with_altloc_all": int(all_.array_length()),
            },
        }

    def gemmi_run() -> dict[str, object]:
        import gemmi

        structure = gemmi.read_structure(str(path))
        residue = structure[0]["A"][0]
        return {
            "default": "keeps every conformer as its own atom with an altloc character",
            "retained": "all conformers; selection is the caller's",
            "result": {
                "atoms_iterated": len(residue),
                "altlocs": sorted({a.altloc for a in residue}),
            },
        }

    rows.append(guarded(case, "molframe", molframe_run))
    rows.append(guarded(case, "biopython", biopython_run))
    rows.append(guarded(case, "biotite", biotite_run))
    rows.append(guarded(case, "gemmi", gemmi_run))
    for each in rows:
        each["probe"] = "synthetic: built from real coordinates, not a corpus entry"
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refined-set", type=Path, required=True)
    arguments = parser.parse_args()
    corpus = json.loads(Path(__file__).with_name("corpus.json").read_text())
    verify_files(arguments.refined_set, [f for e in corpus["entries"] for f in e["files"].values()])
    verify_files(Path(), corpus["bench"])
    table = (
        case_identifiers()
        + case_insertion_codes(arguments.refined_set)
        + case_hydrogens(arguments.refined_set)
        + case_ensemble()
        + case_altloc()
    )
    json.dump({"versions": versions(), "table": table}, sys.stdout, indent=1, default=str)
    print()


if __name__ == "__main__":
    main()
