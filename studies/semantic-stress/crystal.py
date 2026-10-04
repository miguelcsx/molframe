"""Which residues contact another molecule, with and without the crystal around them.

An "interface residue" here is a residue of the deposited unit with an atom that touches
an atom of a different chain instance. Under the asymmetric unit that means another chain
of the file; under ``crystal:<radius>`` it also means a symmetry mate. The audit varies
the system (``assembly``) together with how a contact is defined, and reports how often a
residue's classification changes and which decision does it.

Only the real structures the repository ships are used, each of which has a unit cell.

Run from the repository root:

    PYTHONPATH=python uv run --no-sync python studies/semantic-stress/crystal.py \
        > studies/semantic-stress/crystal.json
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np

import molframe
from molframe import analysis, audit

warnings.simplefilter("ignore")

BENCH = Path("crates/molframe-bench/data")
ENTRIES = ("1crn", "1ubq", "4hhb")
RADIUS = 5.0
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

SPACE = (
    audit.PolicySpace()
    .vary(
        "assembly",
        ["asymmetric_unit", f"crystal:{RADIUS}"],
        rationale="the deposited unit is what the file holds; the crystal is what it was measured in",
        evidence="every entry here has a unit cell and space-group operators",
    )
    .vary(
        "contact_def",
        ["distance:0.0", "distance:0.5"],
        rationale="no slack, and the library's own default of half an angstrom",
        evidence="ContactDefinition::default",
    )
    .vary(
        "vdw_radii",
        ["bondi", "charmm"],
        rationale="two of the four radius sets that ship with the library",
        evidence="molframe-chem radius tables",
    )
)
PLAN = SPACE.plan()
# 4HHB's extent makes its crystal search exceed the default candidate-image ceiling; raising it
# is a stated resource decision, not a change to what is computed.
CONTEXT = molframe.ExecutionContext(image_limit=200_000_000)


def instance_of(origin: np.ndarray, chain_of: np.ndarray) -> np.ndarray:
    """The chain instance each analysed atom belongs to.

    Instances are contiguous blocks of one chain's atoms, in source order, so a new instance
    starts wherever the source index stops increasing by one or the source chain changes.
    """
    sources = chain_of[origin]
    breaks = np.concatenate(
        ([True], (origin[1:] != origin[:-1] + 1) | (sources[1:] != sources[:-1]))
    )
    return np.cumsum(breaks) - 1


def residue_labels(structure: molframe.Structure) -> tuple[list[str], np.ndarray]:
    """Per input atom: a residue label (author chain, number, insertion code) and its chain.

    Only amino-acid residues are labelled: waters, ions and ligands are not part of the
    interface this study asks about, and an unlabelled atom takes no part in it.
    """
    labels = [""] * structure.atom_count
    chains = np.zeros(structure.atom_count, dtype=np.int64)
    for chain_index in range(len(structure.models[0].chains)):
        chain = structure.models[0].chains[chain_index]
        for position in range(len(chain.residues)):
            residue = chain.residues[position]
            if residue.name not in AMINO_ACIDS:
                continue
            label = f"{chain.auth_label}|{residue.auth_number}|{residue.insertion_code or ''}"
            for atom in range(len(residue.atoms)):
                labels[residue.atoms[atom].index] = label
                chains[residue.atoms[atom].index] = chain_index
    return labels, chains


def run_entry(code: str) -> dict:
    structure = molframe.read(BENCH / f"{code}.cif")
    labels, chain_of = residue_labels(structure)
    atoms = structure.atom_count

    def interface(analysed) -> list[str]:
        origin = np.asarray(analysed.atom_origin, dtype=np.int64)
        instance = instance_of(origin, chain_of)
        first = np.asarray(analysed.value["first"], dtype=np.int64)
        second = np.asarray(analysed.value["second"], dtype=np.int64)
        polymer = np.array([bool(label) for label in labels])
        crossing = (
            (instance[first] != instance[second]) & polymer[origin[first]] & polymer[origin[second]]
        )
        # The unit's side of each crossing: the atom of the deposited copy, which comes first.
        unit = np.minimum(first[crossing], second[crossing])
        in_unit = origin[unit][unit < atoms]
        return sorted({labels[i] for i in in_unit})

    result = audit.run(
        PLAN,
        lambda policy: analysis.contacts_by_definition(structure, policy=policy, context=CONTEXT),
        metric="set",
        project=interface,
    )
    sets = [set(interface(run)) for run in result.runs]
    residues = len({label for label in labels if label})  # amino-acid residues only
    by_assembly = {}
    for policy, found in zip(result.policies, sets, strict=True):
        by_assembly.setdefault(policy.assembly, []).append(len(found))
    return {
        "id": code,
        "atoms": atoms,
        "residues": residues,
        "universes": len(sets),
        "interface_residues_by_assembly": {k: sorted(v) for k, v in by_assembly.items()},
        "agreement_with_first": result.agreement_with_first,
        "effects": {effect.field: effect.mean_change for effect in result.effects or []},
        "shapley": {share.name: share.share for share in result.shapley or []},
        "interactions": {f"{i.first}x{i.second}": i.share for i in result.interactions or []},
        "balanced": result.balanced,
    }


def main() -> None:
    rows = []
    for code in ENTRIES:
        rows.append(run_entry(code))
        sys.stderr.write(f"{code}\n")
    json.dump(
        {
            "radius": RADIUS,
            "plan": {
                "universes": PLAN.cost,
                "decisions": [
                    {
                        "field": decision.field,
                        "uncertainty": decision.uncertainty,
                        "rationale": decision.rationale,
                        "evidence": decision.evidence,
                    }
                    for decision in PLAN.decisions
                ],
            },
            "entries": rows,
        },
        sys.stdout,
        indent=1,
    )
    print()


if __name__ == "__main__":
    main()
