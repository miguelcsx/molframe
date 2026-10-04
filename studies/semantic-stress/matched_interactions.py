"""Matched residue-level vdW fingerprints on original deposited structures.

ProLIF and MolFrame receive identical selected atoms, coordinates and radii.
Only VdWContact is evaluated: bond-free RDKit transfer is valid for this geometry
endpoint, and must never be used as a chemical fingerprint or hydrogen-bond test.
The fixed endpoint is standard amino-acid residues contacting each non-water,
non-polymer residue with at least five heavy atoms in the deposited first model.
This is a purposive integration panel, not a population or biological truth test.
"""

import argparse
import hashlib
import json
from importlib import metadata
from pathlib import Path

import numpy as np

import molframe
from molframe import analysis, chemistry

AMINO_ACIDS = frozenset(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)


def rdkit_transfer(structure, indices) -> object:
    """Preserve residue indices and coordinates; no chemical bonds are inferred."""
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    editable = Chem.RWMol()
    coordinates = np.asarray(structure.coordinates)
    conformer = Chem.Conformer(len(indices))
    for number, origin in enumerate(indices):
        atom = structure.atoms[int(origin)]
        transferred = Chem.Atom(atom.atomic_number)
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(atom.residue.name)
        info.SetResidueNumber(atom.residue.index + 1)
        info.SetName(f"{atom.name:>4}")
        transferred.SetMonomerInfo(info)
        transferred.SetUnsignedProp("mapindex", number)
        editable.AddAtom(transferred)
        conformer.SetAtomPosition(number, Point3D(*map(float, coordinates[int(origin)])))
    molecule = editable.GetMol()
    molecule.AddConformer(conformer)
    molecule.UpdatePropertyCache(strict=False)
    return molecule


def entity_kind(structure, residue) -> str | None:
    """Use deposited entity types; modified polymer residues are not ligands."""
    index = residue.chain.entity
    return structure.entities[index].kind if index is not None else None


def prolif_transfer(structure, indices) -> object:
    """Explicit residue grouping avoids ProLIF splitting bond-free atoms into fragments."""
    import prolif
    from prolif.residue import Residue

    groups = {}
    for i in indices:
        groups.setdefault(structure.atoms[int(i)].residue.index, []).append(int(i))
    residues = [Residue(rdkit_transfer(structure, atoms)) for atoms in groups.values()]
    molecule = prolif.Molecule(rdkit_transfer(structure, indices), residues=residues)
    if sum(residue.GetNumAtoms() for residue in molecule.residues.values()) != len(indices):
        msg = "ProLIF transfer lost atoms during residue grouping"
        raise ValueError(msg)
    return molecule


def native_pocket(structure, first, second, protein, ligand) -> set[int]:
    """Map either contact orientation back to original protein residue indices."""
    found = set()
    for a, b in zip(first, second, strict=True):
        if a in protein and b in ligand:
            found.add(structure.atoms[int(a)].residue.index)
        elif b in protein and a in ligand:
            found.add(structure.atoms[int(b)].residue.index)
    return found


def run_configuration(path, hydrogens, radii, slack) -> dict:
    """Compare identical residue-contact questions, retaining differences explicitly."""
    import prolif

    structure = molframe.read(path, options=molframe.ReadOptions(digest_input=True))
    policy = molframe.AnalysisPolicy(
        hydrogens=hydrogens, vdw_radii=radii, contact_def=f"distance:{slack}"
    )
    result = analysis.contacts_by_definition(structure, policy=policy)
    if result.value is None:
        return {"refused": "indeterminate", "policy": [hydrogens, radii, slack]}
    origins = np.asarray(result.atom_origin, dtype=np.int64)
    protein = [
        i
        for i in origins
        if structure.atoms[int(i)].residue.name in AMINO_ACIDS
        and entity_kind(structure, structure.atoms[int(i)].residue) == "polymer"
    ]
    groups = {}
    for i in origins:
        atom = structure.atoms[int(i)]
        residue = atom.residue
        if entity_kind(structure, residue) == "non_polymer":
            groups.setdefault(residue.index, []).append(int(i))
    ligands = {
        key: atoms
        for key, atoms in groups.items()
        if sum(structure.atoms[i].atomic_number > 1 for i in atoms) >= 5
    }
    if not protein or not ligands:
        return {"no_endpoint": True, "policy": [hydrogens, radii, slack]}
    radius_values = np.asarray(chemistry.vdw_radii(structure, radii=radii))
    overrides = {}
    for i in protein + [i for atoms in ligands.values() for i in atoms]:
        if not np.isfinite(radius_values[i]):
            return {"refused": "missing radius", "policy": [hydrogens, radii, slack]}
        overrides[structure.atoms[i].element] = float(radius_values[i])
    fingerprint = prolif.Fingerprint(
        ["VdWContact"], parameters={"VdWContact": {"tolerance": slack, "vdwradii": overrides}}
    )
    receptor = prolif_transfer(structure, protein)
    protein_set = set(map(int, protein))
    first = origins[np.asarray(result.value["first"], dtype=np.int64)]
    second = origins[np.asarray(result.value["second"], dtype=np.int64)]
    comparisons = []
    for residue_index, atoms in sorted(ligands.items()):
        ligand_set = set(atoms)
        native = native_pocket(structure, first, second, protein_set, ligand_set)
        ligand = prolif_transfer(structure, atoms)
        reference = fingerprint.generate(ligand, receptor, residues="all")
        found = {protein_id.number - 1 for (_, protein_id), bits in reference.items() if bits.any()}
        comparisons.append(
            {
                "ligand_residue": residue_index,
                "component": structure.residues[residue_index].name,
                "molframe": sorted(native),
                "prolif": sorted(found),
                "only_molframe": sorted(native - found),
                "only_prolif": sorted(found - native),
            }
        )
    return {"policy": [hydrogens, radii, slack], "comparisons": comparisons}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--ids", nargs="+")
    arguments = parser.parse_args()
    manifest = json.loads(arguments.corpus.read_text())
    root = arguments.corpus.parent
    rows = []
    for record in manifest["files"]:
        if not record["path"].endswith(".cif") or "error" in record:
            continue
        path = root / record["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            msg = f"frozen input changed: {path}"
            raise ValueError(msg)
        if arguments.ids and path.stem not in arguments.ids:
            continue
        for hydrogens in ("explicit_only", "exclude"):
            for radii in ("bondi", "charmm", "alvarez", "amber_united"):
                if hydrogens == "explicit_only" and radii == "amber_united":
                    continue
                for slack in (0.0, 0.25, 0.5):
                    try:
                        row = run_configuration(path, hydrogens, radii, slack)
                    except (ValueError, RuntimeError) as error:
                        row = {"policy": [hydrogens, radii, slack], "refused": str(error)}
                    rows.append({"id": path.stem, **row})
    print(
        json.dumps(
            {
                "molframe": molframe.__version__,
                "corpus_sha256": hashlib.sha256(arguments.corpus.read_bytes()).hexdigest(),
                "versions": {name: metadata.version(name) for name in ("prolif", "rdkit", "boba")},
                "scope": "Matched vdW residue fingerprints, shared selection/radii/coordinates; no chemical or biological ground truth claim",
                "runs": rows,
            },
            indent=1,
        )
    )


if __name__ == "__main__":
    main()
