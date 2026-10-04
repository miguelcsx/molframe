# Semantic stress: how far does an answer move when a defensible decision changes?

This directory holds a small, frozen, reproducible study. It asks one question of
real structures: *for an analysis whose inputs are fixed, how much of the answer
depends on decisions that the file does not make and the analyst usually does not
notice?* It is a review of how MolFrame's policy and audit layer behaves on real
data, and it records where that layer fell short while the study was being built.
It is not a benchmark of other tools and it does not say which choice is correct.

The original tables below are read from the historical JSON files beside this
README. The 2026-10-04 review at the end records a fresh run with a corrected
decision space; the historical 24-universe pocket tables do not describe that run.

## What was run

| Script | Output | What it does |
|---|---|---|
| `build_corpus.py` | `corpus.json` | Freezes the entries, strata and SHA-256 of every file. |
| `tool_semantics.py` | `tool_semantics.json` | The same structural condition read by MolFrame, Biopython, Biotite and Gemmi. |
| `mdanalysis_semantics.py` | `mdanalysis_semantics.json` | The same, for MDAnalysis, which reads PDB but not mmCIF. |
| `pockets.py` | `pockets.json` | Ligand pockets of 60 complexes under 24 universes of three decisions. |
| `burial.py` | `burial.json` | How buried the ligands of the same 60 complexes are under radii and hydrogens, with a forbidden pair. |
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
been assembled. It needs a prespecified sampling design, identifiers, archived
versions and validation metadata before a population study can be run.

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

## How buried a ligand is (`burial.json`)

A ligand atom is *buried* when its solvent-accessible area in the complex is under
1 Å², and a ligand is *mostly buried* when more than half of its atoms are. The threshold
and the fraction are parameters of the question, fixed here by me and not varied; a
different pair would give different numbers. The decisions varied are the radii the
atoms are given (`vdw_radii`: bondi, charmm, alvarez, amber_united) and whether the
file's hydrogens occlude (`hydrogens`: explicit_only, exclude).

United-atom radii already contain their hydrogens, so `amber_united` with modelled
hydrogens counts them twice. That pair is forbidden in the space with that reason, which
makes the plan 7 universes of the 8, not the whole product, and the attribution below is
Shapley's rather than the additive split's.

Atoms whose element has no radius in a set (a zinc or calcium ion under Bondi) have no
area and do not occlude; the analysis is *partial* under the default `missing_atoms`
rule and says which atoms. 59 of 60 complexes ran (1C5S is the same occupancy refusal as
above). One complex has a ligand atom without a radius in some universes.

| Quantity (59 complexes, 7 universes each) | Value |
|---|---|
| Fraction of ligand atoms buried, mean over complexes of the midpoint of its range | 0.80 |
| Range of that fraction across universes, median (max) | 0.09 (0.32) |
| Complexes where "mostly buried" is true in some universes and false in others | **6 of 59 (10 %)**: 2DRC, 2W5G, 3CD5, 3ZK6, 4AVS, 4OVG |
| Ligand exposed area, relative change from one decision, median: hydrogens / radii | 6 % / 10 % |
| Shapley share of the exposed-area variation: radii / hydrogens | 0.68 / 0.28 |
| Shapley share of the buried-fraction variation: hydrogens / radii | 0.52 / 0.36 |

For most ligands the verdict "mostly buried" does not depend on these choices, because
most ligands are far from the line. The six that flip are near it, and for them the
sentence "the ligand is buried" is a statement about the radius table as much as about the
complex. Which table is closer to the truth is not something this audit can say. The
attribution depends on what is asked: the area as a number moves mainly with the radii,
while the count of buried atoms against a threshold moves more with the hydrogens.

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
- Hydrogen bonds, salt bridges, contact maps and structural comparison have each been
  verified separately against independent calculations in the test suite but have not been
  swept at scale here. Protein-ligand hydrogen bonds need a chemical dictionary for each
  ligand, which this corpus does not ship.

## Reproducing

From the repository root, build the current extension in the development shell:

```bash
nix develop --command sh -c \
    'unset UV_PYTHON; VIRTUAL_ENV="$PWD/.venv" maturin develop --release --uv'
```

Then run the study in that shell:

```bash
nix develop --command sh -c 'PYTHONPATH=python uv run --no-sync \
    --with biopython==1.88 --with biotite==1.7.1 --with gemmi==0.7.5 \
    python studies/semantic-stress/pockets.py \
        --refined-set ../ligare/data/refined-set > /tmp/molframe-pockets.json'
```

The refined set is not redistributed here. `corpus.json` records the SHA-256 of every file
used, so a copy that differs is detected rather than silently used.

## Review and rerun, 2026-10-04

The runners now verify SHA-256 against the frozen manifest before analysis.
The three governed studies request input digests, and their outputs record the
corpus digest and software version. A digest of the combined PDB identifies the
derived input; the frozen manifest identifies its source protein and ligand.
Combining SDF coordinates into PDB rounds them to 0.001 Å, so these studies
measure the derived representation, not the original SDF precision.

The pocket plan now excludes explicit hydrogens with united-atom radii, as the
burial plan already did. It has **21 universes, three excluded**, and is
constrained. Its output withholds additive interaction and higher-order shares;
the old 16% interaction statement above is historical and must not be cited for
the corrected plan.

| Fresh check | Result |
|---|---|
| Ligare refined-set protein reads | 5,316 / 5,316 |
| Ligare refined-set ligand SDF reads, after the parser fix | 5,316 / 5,316 |
| Frozen 60-entry source hashes | All match |
| Pockets, 21 universes per analysable complex | 59 / 59 change; one occupancy refusal (1C5S) |
| Pocket pairwise Jaccard distance, mean over complexes | 0.266716 |
| Distinct pockets, median (range) | 11 (4–17) |
| Pocket Shapley shares, mean: contact definition / radii / hydrogens | 0.446812 / 0.352927 / 0.200261 |
| Independent NumPy pocket geometry | 126 / 126 agree: six complexes × 21 universes |
| Burial, seven universes | 6 / 59 flip; one occupancy refusal |
| Crystal system-definition probe | Reproduces 1CRN, 1UBQ and 4HHB ranges above |
| Biopython / Biotite / Gemmi and MDAnalysis semantics probes | Rerun at the recorded versions |

The independent pocket check uses NumPy distances and the same MolFrame radius
tables; it validates geometry, filtering and origin mapping, not independent
chemical-radius correctness. Input reads alone do not validate all analyses on
all 5,316 complexes. The multiverse studies still use the 60-entry pilot.

Operationally the pocket projection includes every receptor-file residue except
`HOH`, including ions and cofactors. The earlier shorthand “protein residues”
must not be used to claim a protein-only endpoint; that would require a separately
validated residue-class mask. The same projection is used for the reference tools.

Reading the full refined set exposed 264 failed SDF reads. V2000 counts and bond
indices are fixed-width fields that become adjacent at three digits; the old
parser split them by whitespace. It now reads counts, coordinates and bonds
by column, with tests for three-digit fields and touching full-width coordinates.
Existing Rust fixtures also used line-continuation escapes that stripped required
leading blanks; they now preserve the actual fixed columns.

Burial flips remain 2DRC, 2W5G, 3CD5, 3ZK6, 4AVS and 4OVG. One different
complex, 1C5Q, lacks ligand radii in some universes. Excluding it gives a
full-ligand-coverage subset of 58, with six flips; report that denominator as
well as the covered-atom analysis. Neither rate estimates prevalence in the PDB.

The local fresh outputs are `/tmp/molframe-{pockets,burial,crystal}-final.json`,
`/tmp/molframe-ligare-check-final.json`, `/tmp/molframe-tool-semantics-current.json`
and `/tmp/molframe-mdanalysis-current.json`. These are temporary run artifacts,
not additional committed generated data. Reproduce and archive the output of a
prespecified study before using any number in a manuscript.

## Matched domain references and original data

The newer comparison asks the same scientific questions as its references:
`matched_interactions.py` compares **residue-level vdW fingerprints with ProLIF**;
`matched_surface.py` compares **ligand SASA and burial with FreeSASA through
RDKit**, using Lee-Richards with 20 slices. The earlier generic library tests
remain reader/geometry checks and do not establish interaction-method superiority.
PLIP and Arpeggio are relevant chemical-fingerprint comparators but were not run;
the point-cloud transfer here is valid only for geometry, not chemical perception.

Inputs and outputs live in `../ligare/data/molframe-research`, outside git.
`acquisition.json` freezes URLs, bytes and SHA-256 for 12 deposited mmCIF files,
12 wwPDB validation reports and the BioLiP-like Q-BioLiP annotation download
(602,181 rows). Seventeen annotation rows refer to panel entries. These curated
binding residues use a geometric criterion and quaternary structures; reconcile
assembly/copy and identifiers before scoring, and do not call them independent
experimental contact truth. Validation XML is downloaded for external QC;
quality-prediction experiments remain pending. The panel is purposive and
contains related proteins; it cannot estimate archive-wide prevalence.

| Matched question | Fresh result |
|---|---|
| ProLIF 2.2.2 / MolFrame vdW residue fingerprints, RDKit 2026.3.6 | 150 ligand fingerprints agree exactly |
| Requested runs | 252: 12 structures × 21 admissible policies |
| No endpoint | 84 runs across four controls (1UBQ, 1CRN, 1BNA, 2LZM) |
| Refused | 21 runs: 4DFR occupancy 1.02; 15 runs: missing ligand radii in 4HHB |
| FreeSASA / MolFrame per-atom Lee-Richards SASA | Four complexes; maximum absolute atom error 1.16e-12 Å² |
| Ligand burial fraction, threshold <1 Å² per atom | Identical in all four FreeSASA comparisons |
| Boba 1.1.2, same calculation for 1HVR | 21 generated/executed universes; all output contents match the direct runner |

The 132 analysable fingerprint runs supply 150 comparisons because some
structures have several ligands. Protein atoms are standard amino-acid residues
in polymer entities; ligands are non-polymer residues with at least five heavy
atoms. Modified polymer residues are not ligands. ProLIF receives explicit
residue groups and verifies atom conservation; its default connected-fragment
splitting is unsuitable for the bond-free geometry transfer. Both references
receive matching selected atoms, float32 coordinates, supplied radii and probe
or tolerance. This isolates kernels, not independent chemistry or selection.
Surface per-atom roundoff tolerance is 1e-8 Å²; the SASA references agree much
more closely. Missing ligand radius coverage is refused rather than silently
changing the ligand denominator. Neither software agreement nor stability is
biological correctness.

Boba's successful run shows that constrained enumeration is already available
in a generic multiverse tool. Its older runner labels any stderr text “error”,
including ProLIF's upstream MDAnalysis deprecation warning. All 21 result files
were independently verified; no performance comparison was inferred.

To rerun the domain comparisons, use the Nix development environment and an
isolated Python environment containing the current MolFrame extension,
`prolif==2.2.2`, `rdkit==2026.3.6`, and `boba==1.1.2`:

```bash
python studies/semantic-stress/matched_interactions.py --corpus ../ligare/data/molframe-research/acquisition.json > matched-interactions.json
python studies/semantic-stress/matched_surface.py --corpus ../ligare/data/molframe-research/acquisition.json > matched-surface.json
```

The manifest is a downloaded snapshot, not bundled test data. Its acquisition
URLs can be fetched again, but changed bytes require an explicit new manifest.
Original-to-analysis identity remains in original residue indices and frozen
input bytes. The new governed audit also refuses mismatching declared estimands
before attribution; this checks declarations, not full biological equivalence.

For the Boba comparison, save this template outside the checkout. Use an absolute
`MOLFRAME_CORPUS` path and include `studies/semantic-stress` in `PYTHONPATH` for
the generated programs. Compile with `boba compile -s template.py --out output`,
then run `boba run --all --dir output/multiverse` from the same directory. Compare
all 21 `code/result-*.json` contents with the direct runner's 1HVR rows, matching
by policy rather than file order.

```python
# --- (BOBA_CONFIG)
{"decisions": [{"var": "hydrogens", "options": ["explicit_only", "exclude"]}, {"var": "radii", "options": ["bondi", "charmm", "alvarez", "amber_united"]}, {"var": "slack", "options": [0.0, 0.25, 0.5]}], "constraints": [{"variable": "slack", "index": 0, "condition": "hydrogens == exclude or radii != amber_united"}, {"variable": "slack", "index": 1, "condition": "hydrogens == exclude or radii != amber_united"}, {"variable": "slack", "index": 2, "condition": "hydrogens == exclude or radii != amber_united"}]}
# --- (END)
import json
import os
from pathlib import Path
from matched_interactions import run_configuration

result = run_configuration(
    Path(os.environ["MOLFRAME_CORPUS"]) / "deposited/1hvr.cif",
    "{{hydrogens}}", "{{radii}}", {{slack}},
)
Path("result-{{_n}}.json").write_text(json.dumps(result))
```

The direct runner verifies the frozen corpus hashes before analysis; verify
that same manifest before the Boba execution. Boba is a runner reference here,
not an independent interaction kernel.


## Experimental interface pilot

`prepare_skempi.py` freezes a deterministic panel from the original SKEMPI 2.0
publisher XLS and a commit-pinned CSV mirror. Source URLs/hashes are retained in
`../ligare/data/molframe-research/experimental/source-{publisher,mirror}.json`.
The mirror's identifiers, affinities, holdout fields, reference and temperature
must agree with the XLS. The XLS has space-split fields; discordance is an explicit
exclusion, never a schema guess used silently. Numeric single-alanine rows are
eligible; censored/nonpositive affinities, unspecified temperatures, WT identity
conflicts and explicitly mismatched WT/mutant temperatures are excluded.

The fourteen downloaded original structures were selected by SHA-256 order,
minimum ten distinct eligible mutations and at most one per connected component
of supplied holdout links/shared protein names. The protocol hash is checked on
every analysis. A sequence screen later tests cross-group similarity; it is not
exhaustive family adjudication. All data/generated reports remain outside git.

`experimental_interfaces.py` matches measured author-chain/residue/insertion/WT
identities and measures twelve heavy-atom vdW contact-degree policies. The fixed
baseline is Bondi + 0.5 Å; the aggregate is their uniform median. Experimental
ddG uses the raw Kd ratio and numeric temperature. Replicate/condition records are
retained and a measured residue is one statistical unit. A per-residue median
summarizes reported assays rather than estimating a common-condition free energy.
Seven inconsistent hotspot labels remain ambiguous, not forced to zero or one.
Thresholds are calibrated outside the entire test metadata group. Three systems
with a single class retain their records but have unavailable AUROC/accuracy.

The corrected panel has 631 measures of 545 residues, 168 completed configurations
and no structure-level refusals. Baseline versus aggregate mean AUROC is
0.7223/0.7197, and held-out balanced accuracy is 0.6755/0.6871 across eleven
systems. Conditional paired bootstrap results retain their limitations: predictions
are not refitted within resamples. These do not establish confirmatory significance
or consistent superiority of aggregation. Some fixed policies outperform it.
204 residue degrees vary; a wrong label-namespace join changes/loses 371 residues
across eight systems. A checked generic runner can implement the same safeguards.

After preliminary results, the assay-condition review removed 32 mixed-temperature
rows and one contradictory experimental WT identity. The same structures remained
selected. Protocol/results before this correction are archived with
`.pre-condition-check.json` suffixes. This correction is disclosed rather than
presented as external preregistration.

`matched_experimental.py` verifies author joins with Gemmi, screens sequence
relatedness using Biopython, and checks interface degree against ProLIF on 1BRS,
1GUA and 4PWX under all twelve policies. Shared atoms/radii isolate geometry.
Chemical perception, experimental construct/state compatibility and a matched
energy-predictor benchmark remain outside this workflow.

Use an isolated environment containing the current MolFrame extension and pinned
references; the XLS extraction uses `xlrd==2.0.2`, the reference workflow uses
Gemmi 0.7.5, ProLIF 2.2.2, RDKit 2026.3.6 and the recorded Biopython version.
SciPy computes rank correlation; NumPy computes descriptive paired bootstraps.
From the Nix development environment:

```bash
python studies/semantic-stress/prepare_skempi.py --inputs ../ligare/data/molframe-research/experimental > protocol.reproduced.json
# Compare with frozen protocol.json; changing inputs/protocol requires a new explicit acquisition manifest.
python studies/semantic-stress/experimental_interfaces.py --protocol ../ligare/data/molframe-research/experimental/protocol.json > results.reproduced.json
python studies/semantic-stress/matched_experimental.py --results ../ligare/data/molframe-research/experimental/results.json > reference-results.reproduced.json
```

The publisher's supplementary-data link may require a renewed signed URL. Preserve
the source bytes, mirror commit, exclusions and manifest when reproducing rather
than silently substituting a changed download. Full numerical outputs and a static
Matplotlib figure (`policy-evidence.png`/`.pdf`) are archived beside the corpus.


`matched_energy.py` adds a descriptive binding-change predictor comparison with
EvoEF 1.1, source commit `6bce56dcbf95e34dfe890bfccdffd0df252b830c`, compiled
outside MolFrame with Nix Clang 21.1.8 (`-O3 -ffast-math`). It uses the same
three reference complexes, author targets and heavy-atom protein input selection.
Repair uses three passes and mutant optimization ten, with an explicit partner
split. Insertion codes are refused by this transfer. Full residue maps are checked
after repair and after each of 39 mutations; model/input/log/binary/library hashes
are retained. Model construction and chemical energetics are specific to EvoEF.

1BRS energy/contact-degree AUROC is 0.7875/0.8000; Spearman is 0.6227/0.6117.
1GUA and 4PWX have only one hotspot class, so AUROC is unavailable. This is a
small descriptive comparison, with unresolved overlap between pilot data and
EvoEF's SKEMPI-derived training. No independent-generalization claim follows.
Energy differences and contact counts have different units; compare their ranks,
not their absolute errors against one another.

```bash
python studies/semantic-stress/matched_energy.py --reference ../ligare/data/molframe-research/experimental/EvoEF-reference/EvoEF --results ../ligare/data/molframe-research/experimental/results.json --output ../ligare/data/molframe-research/experimental/energy > energy-results.reproduced.json
```

The external reference checkout has the original license and source. MolFrame
contains only a hand-written CLI adapter, with no copied reference implementation.


### Rust/Python corpus parity (0.5.0)

The optional native binding test
`python_and_rust_read_the_same_external_ligand_corpus` compared all 5,316 SDFs in
Ligare's refined-set. Atom counts and all coordinates agreed exactly. To rerun,
set `MOLFRAME_PARITY_CORPUS` to the acquired refined-set, `PYO3_PYTHON` to the
project Python interpreter and `PYTHONPATH` to its site-packages, then run inside
`nix develop`: `cargo test -p molframe-py python_and_rust_read_the_same_external_ligand_corpus -- --ignored --nocapture`.
This checks transfer/reading parity; independent scientific references remain
the separate experiments above.
