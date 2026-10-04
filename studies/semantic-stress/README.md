# Semantic stress: how far does an answer move when a defensible decision changes?

This directory holds a small, frozen, reproducible study. It asks one question of
real structures: *for an analysis whose inputs are fixed, how much of the answer
depends on decisions that the file does not make and the analyst usually does not
notice?* It is a review of how MolFrame's policy and audit layer behaves on real
data, and it records where that layer fell short while the study was being built.
It is not a benchmark of other tools and it does not say which choice is correct.

Everything below is read from the JSON files beside it. Nothing is recalled.

## What was run

| Script | Output | What it does |
|---|---|---|
| `build_corpus.py` | `corpus.json` | Freezes the entries, strata and SHA-256 of every file. |
| `tool_semantics.py` | `tool_semantics.json` | The same structural condition read by MolFrame, Biopython, Biotite and Gemmi. |
| `mdanalysis_semantics.py` | `mdanalysis_semantics.json` | The same, for MDAnalysis, which reads PDB but not mmCIF. |
| `pockets.py` | `pockets.json` | Ligand pockets of 60 complexes under 24 universes of three decisions. |
| `crystal.py` | `crystal.json` | Interface residues of three crystals with and without the crystal around them. |

Library versions are recorded in the outputs: Biopython 1.88, Biotite 1.7.1,
Gemmi 0.7.5, MDAnalysis 2.10.0.

## The corpus, and what it cannot show

- **60 protein–ligand complexes** from the PDBbind v2020 refined set, twelve from each
  of five strata (`control`, `insertion-codes`, `multi-chain`, `resolution-high`,
  `resolution-mid`), chosen by ascending `sha256(stratum:id)`. No random seed and no
  directory-order dependence. The 5,316-entry set is summarised in `corpus.json`.
- **Three crystals** that ship with the repository (1CRN, 1UBQ, 4HHB), plus 1AON, the
  32-model NMR entry 2M7C, and entries 1C1U and 6D9X for the tool table.

**This is not the corpus the question needs.** PDBbind files are prepared:
none of the 5,316 protein files has a `CRYST1` record (checked), so there is no unit
cell and no assembly; none of the 60 chosen entries has an alternate location
(`altloc_atoms` is 0 in `corpus.json`); and hydrogens have been added by a program. So the
strata that most need testing (alternate locations, several assemblies, missing
atoms and residues, ensembles, author/label divergence, symmetry contacts in
a dataset with real variety) are absent from the 60 complexes. Altlocs appear only
as a *synthetic* two-conformer probe built from 1UBQ's real coordinates, labelled
as such in the table. A corpus of original mmCIF entries for those strata has not
been downloaded; it needs the person who owns this machine to approve a list of
identifiers, URLs and sizes, which this study does not yet have.

## Ligand pockets (`pockets.json`)

For each complex the pocket is the set of protein residues with an atom in contact
with a ligand atom. It is recomputed under every combination of:

| Decision | Class | Values |
|---|---|---|
| `hydrogens` | interpretive | the file's hydrogens count / excluded |
| `vdw_radii` | algorithmic | bondi, charmm, alvarez, amber_united |
| `contact_def` | algorithmic | distance slack of 0.0, 0.25, 0.5 Å over the sum of radii |

24 universes per complex. 59 of 60 complexes ran; **1C5S was refused (`E3009`,
occupancy outside zero to one)**, because the file carries an occupancy of 1.35 and
the executor does not guess what that means. The three other tools read it without
comment. The refusal is counted as a result and is excluded from every figure below.

| Quantity (59 complexes) | Value |
|---|---|
| Pocket differs between at least two universes | **59 / 59** |
| Distinct pockets per complex (median, range) | 12 (4 – 20) of 24 universes |
| Pocket size, union over universes (median) | 16 residues |
| Residues in the union that are not in every universe (mean fraction) | 0.68 |
| Mean pairwise Jaccard distance between universes | 0.26 (max 0.95) |
| Agreement of a universe with the first one | 0.10 |

Changing a *single* decision, the rest held:

| Decision changed | Mean Jaccard change, median over complexes | Range | Complexes with change > 0.1 |
|---|---|---|---|
| `hydrogens` | 0.17 | 0.03 – 0.38 | 78 % |
| `vdw_radii` | 0.19 | 0.04 – 0.31 | 92 % |
| `contact_def` | 0.26 | 0.05 – 0.42 | 95 % |

Shapley attribution of the variation (mean over complexes): `contact_def` 0.46,
`vdw_radii` 0.35, `hydrogens` 0.19. By class, **algorithmic 0.80, interpretive
0.20**. Mean interaction shares: `vdw_radii × contact_def` 0.16 (the two together
matter more than either alone), `hydrogens × contact_def` 0.07, `hydrogens × vdw_radii`
0.06, higher order 0.05. The ordering by stratum is flat (mean pairwise distance
0.24 – 0.28), so no stratum is a more fragile place to ask the question.

### What this does and does not mean

- "100 % of pockets change" is a statement about **this universe set**, whose members
  are not equally plausible. Four radius tables and three slacks include choices a
  group would not make together. It says the pocket is a function of these choices, not
  that any published pocket is wrong.
- The audit measures sensitivity. It cannot say which universe is closer to the
  truth, and the study does not claim one is.
- The spread is mostly the size of the *question* (a contact is a threshold on
  distance, and a pocket is the residues past it), not unusual to this library. The
  value is that the fraction is now measured per decision, and that the interaction of
  radii and slack, which none of the two alone shows, is 16 % of the variation.
- 60 complexes, one dataset, one set of three decisions: the table above is a
  demonstration that the machinery runs on real data, not an estimate for the PDB.

### Where other tools' defaults fall

Each tool was asked for the protein residues within 4 Å of any ligand atom, all
atoms, no hydrogens filter. The three tools returned the **same** pockets (they ask
one geometric question), so one row stands for all.

| | Value |
|---|---|
| The tool's pocket equals one of MolFrame's 24 universes | 6.8 % of complexes |
| Mean Jaccard overlap with the *nearest* universe | 0.88 |
| Mean Jaccard overlap with the first universe | 0.52 |

A plain distance cutoff and a van der Waals contact are different definitions. The
result is that a "pocket" from a cutoff is rarely one a radius-based definition
reproduces exactly, even though it is usually near one.

## Interface residues of crystals (`crystal.json`)

An interface residue is an amino-acid residue of the deposited unit with an atom that
touches an atom of a different chain instance. The audit varies the system together
with how a contact is defined: `assembly` (asymmetric unit, or `crystal:5.0`),
`contact_def` (0.0 or 0.5 Å) and `vdw_radii` (bondi or charmm), 8 universes each.

| Entry | Amino-acid residues | Asymmetric unit | Crystal (5 Å) | Shapley share: `assembly` |
|---|---|---|---|---|
| 1CRN (one chain) | 46 | 0 | 12 – 31 | 0.91 |
| 1UBQ (one chain) | 76 | 0 | 12 – 45 | 0.87 |
| 4HHB (four chains) | 574 | 22 – 69 | 44 – 142 | 0.43 |

For a single-chain entry the asymmetric unit has no interface at all, and the
crystal has an interface of a quarter to a half of the residues; the system, not the
contact rule, is nearly all of the variation. For 4HHB, whose file already holds the
tetramer, the crystal roughly doubles the interface and the contact rule matters
about as much as the system.

Crystal contacts are not biological interfaces. 1CRN is a monomer in solution.
The point of the row is how much a "surface" or "interface" analysis would change if
it silently used the asymmetric unit, not that crystal contacts are interactions.

## What each tool does with the same input (`tool_semantics.json`, `mdanalysis_semantics.json`)

Condition → default behaviour → what the tool retains → what came out. A default is a
choice, not an error; the table is here because the choices differ and are silent.

| Input condition | Tool | Default behaviour | Information retained | Result |
|---|---|---|---|---|
| `auth_seq_id` ≠ `label_seq_id` (1AON) | MolFrame | selectors read author ids | both, as columns | 5 atoms either way |
| | Biopython | `id[1]` is the author id; chain is `auth_asym_id` | label ids only in the raw `MMCIF2Dict` | 5 by author, **27 by label**: the same number names a different residue |
| | Biotite | `use_author_fields=True` | either, by re-reading; not both in one array | 5 / 5 |
| | Gemmi | `seqid` is the author id | label id kept per residue | 5 / 5 |
| | MDAnalysis | no mmCIF reader | n/a | n/a |
| Insertion codes (1C1U, 574 residues) | MolFrame | the residue is a topology row; the code is a column | insertion code | 574 |
| | Biopython | residue id is (hetero flag, number, code) | code | 574 |
| | Biotite | `res_id` is the number; code is an annotation | code | 574, but 531 distinct (chain, number) |
| | Gemmi | `ResidueId` is (number, code) | code | 574 |
| | MDAnalysis | run of atoms with equal (resid, resname, segid) | `icodes` attribute | 574 residues, 286 distinct `resid`, 330 with code |
| 32-model NMR ensemble (2M7C) | MolFrame | the first model unless the policy names another or `all` | all 32 | 300 atoms/model |
| | Biopython | iteration yields all models | all 32 | 300 |
| | Biotite | `model=None` gives a stack | all 32 | 300 |
| | Gemmi | all models | all 32 | 300 |
| | MDAnalysis | models become frames, first on the Universe | all 32 | 32 frames |
| Alternate conformations (**synthetic** probe from 1UBQ MET1) | MolFrame | one self-consistent label per region, recorded | all conformers, with label and occupancy | 6 of 8 atoms analysed by default, status `ambiguous` |
| | Biopython | highest-occupancy conformer in `DisorderedAtom` | all, inside it | 6 atoms iterated |
| | Biotite | `altloc="first"` | all with `altloc="all"` | 6 by default, 8 with `all` |
| | Gemmi | every conformer is its own atom | all | 8 |
| | MDAnalysis | every conformer is its own atom | `altLocs`, occupancy | 8 |
| Explicit hydrogens (6D9X): heavy-atom SASA | MolFrame | hydrogens are atoms that occlude | all atoms | 9,687 with H, 12,268 without Å² |
| | Biotite | radius table (ProtOr) in which a hydrogen contributes nothing | all atoms | 12,222 either way |
| | Biopython, Gemmi | no comparable default | | not run |
| | MDAnalysis | hydrogens are atoms; `not name H*` removes by *name* | all atoms | 4,901 atoms; 2,289 hydrogens by element, **1,754 by name prefix** |

Observations worth keeping from building this table:

- MolFrame's heavy-atom SASA without hydrogens (12,268) agrees with Biotite's (12,222)
  to 0.4 %; with the file's hydrogens it is 21 % lower. A user who runs both tools
  without removing hydrogens compares two different quantities.
- In MDAnalysis a name-based hydrogen filter misses a quarter of the hydrogens of
  this file (1,754 of 2,289), because many hydrogens are named with a leading digit (`1HB`), a convention these files use.
- Under MolFrame's default, a *selection* of the altloc probe still returns all 8
  atoms; only governed analyses resolve conformers. Whether that split is
  documented enough is a separate question; a user who selects and then measures gets a different answer from one who analyses
  directly. It is the sharpest edge found, and a candidate for a policy-aware selection.

## What this review changed in MolFrame

The study was run while the audit layer was being finished, and each of these was a defect it exposed:

- Results were indexed in the *analysed* system, not in the input, so a projection
  could not be tied to a residue. `Analysis::atom_origin` now maps each analysed atom
  to its input atom (also in Python).
- A crystal view omitted the inverse placement of a mutual contact: 7 versus 10 contacts in
  a brute-force NumPy cross-check. Fixed, with the check as a test.
- 4HHB's crystal search exceeded the default candidate-image ceiling. The ceiling
  moved to `ExecutionContext` (`image_limit`) so that raising it is a stated resource
  decision.
- The audit reported a constrained plan as balanced. The plan now owns that fact.
- A vdW-radius or contact-definition sweep returned a perfectly stable answer from an
  analysis that never read either. `contacts_by_definition` now reads both and the
  audit refuses a sweep over a decision nothing read (`E6103`).

## Limits, stated plainly

- No original mmCIF with real altlocs, several assemblies, missing atoms or
  author/label divergence beyond 1AON is in the corpus yet.
- The universes are chosen, not sampled; their plausibility is not weighted.
- `Fragility` is not related to model quality, resolution or B-factor here. That needs
  more entries than 60 complexes of a prepared set.
- Interfaces, hydrogen bonds, salt bridges, SASA, contact maps and structural comparison
  have each been verified separately against independent calculations in the test
  suite. They have not yet been swept at scale in this study.

## Reproducing

From the repository root, with the extension built (`maturin develop --release`):

```bash
PYTHONPATH=python uv run --no-sync \
    --with biopython==1.88 --with biotite==1.7.1 --with gemmi==0.7.5 \
    python studies/semantic-stress/pockets.py \
        --refined-set /path/to/refined-set > studies/semantic-stress/pockets.json
```

The refined set is not redistributed here. `corpus.json` records the SHA-256 of every file
used, so a copy that differs is detected rather than silently used.
