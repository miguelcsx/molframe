"""Cut the bundled core component set from the full wwPDB dictionary.

usage: extract_core.py components.cif.gz snfg.rs out.cif release-date

Each kept block retains its identity, atoms (without coordinates) and bonds.
"""  # noqa: INP001

import gzip
import re
import sys
from pathlib import Path

CORE = """
ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
ASX GLX UNK SEC PYL MSE SEP TPO PTR CSO HYP MLY M3L KCX LLP CME CSD PCA FME ACE
NH2 NME ALY OCS CGU MEN HIC SAC TYS CSS CSX SMC DAL DAR DSG DAS DCY DGN DGL DHI
DIL DLE DLY MED DPN DPR DSN DTH DTR DTY DVA ORN ABA NLE CAS
A C G U I DA DC DG DT DU DI N DN PSU 5MC OMC OMG 1MA 2MG 5MU H2U M2G 7MG OMU
1MG 5BU 5CM 6MA 8OG A2M 3DR PUR 70U UR3 5FU 4SU
HOH DOD
NA K CL CA MG ZN FE FE2 MN CU CU1 CO NI CD HG BR IOD F LI RB CS SR BA AL AU PT
TL MO W PB AG YB GD SM EU LA CE NH4
SO4 PO4 NO3 ACT GOL EDO PEG PGE MPD DMS BME TRS MES EPE FMT IPA CIT TAR SCN AZI
IMD CO3 OH PGO PG4 1PE TLA DTT MLI MLA OXL ACY BCT CAC CXS CYN
HEM HEC HEA HEB FAD FMN NAD NAP NDP NAI ATP ADP AMP ANP GTP GDP GNP GSP UTP
UDP CTP CDP TTP SAM SAH COA ACO PLP PMP TPP BTN FES SF4 F3S CLA RET H4B BH4
THF FOL NMN FNR BOG OLC LDA SDS DAO MYR PLM STE OLA
""".split()  # noqa: SIM905
DROP = set("15P CNC COB B12 BCL PHO LMT HEO FDA MGD BCR CHD MTX PQQ PE4 P6G".split())  # noqa: SIM905

TOKEN = re.compile(r"""'(?:[^']|'(?!\s))*'|"(?:[^"]|"(?!\s))*"|\S+""")
META = ["id", "name", "type", "formula", "mon_nstd_parent_comp_id", "one_letter_code"]
ATOM = [
    "atom_id",
    "alt_atom_id",
    "type_symbol",
    "charge",
    "pdbx_aromatic_flag",
    "pdbx_leaving_atom_flag",
    "pdbx_stereo_config",
]
BOND = ["atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag", "pdbx_stereo_config"]
KEPT = {"chem_comp_atom": ATOM, "chem_comp_bond": BOND}


class Block:
    """Incremental reader of one dictionary block."""

    def __init__(self, header: str) -> None:
        """Start a block whose header line is `header`."""
        self.out = [header]
        self.meta: dict[str, str] = {}
        self.items: list[str] = []
        self.rows: list[list[str]] = []
        self.single: dict[str, str] = {}
        self.category = ""
        self.state = ""
        self.last = ""
        self.text: list[str] | None = None

    def flush(self) -> None:
        """Emit the pending category if it is one that is kept."""
        if self.single:
            self.items = list(self.single)
            self.rows = [list(self.single.values())]
            self.single = {}
        wanted = KEPT.get(self.category)
        if wanted and self.rows:
            names = [name for name in wanted if name in self.items]
            self.out.append("loop_\n")
            self.out.extend(f"_{self.category}.{name}\n" for name in names)
            columns = [self.items.index(name) for name in names]
            self.out.extend(" ".join(row[c] for c in columns) + "\n" for row in self.rows)
        self.items, self.rows, self.category = [], [], ""

    def meta_text(self, line: str) -> None:
        """Read a metadata value that continues on following lines."""
        if line.startswith(";"):
            if self.text is None:
                self.text = [line[1:].rstrip("\n")]
            else:
                self.meta[self.last] = '"' + " ".join(self.text).replace('"', "'") + '"'
                self.text = None
        elif self.text is not None:
            self.text.append(line.rstrip("\n"))
        else:
            self.meta[self.last] = line.strip()

    def single_item(self, line: str) -> None:
        """Read one `_category.item value` line of a single-row category."""
        category, _, item = line.strip()[1:].partition(".")
        item, _, value = item.partition(" ")
        if self.state != "single" or self.category != category:
            self.flush()
            self.category, self.state = category, "single"
        self.single[item] = value.strip()

    def loop_row(self, line: str) -> None:
        """Read one data row of a loop."""
        self.state = "data"
        if self.category not in KEPT:
            return
        tokens = TOKEN.findall(line)
        if len(tokens) != len(self.items):
            message = f"unsupported row layout: {line!r}"
            raise SystemExit(message)
        self.rows.append(tokens)

    def feed(self, line: str) -> None:
        """Consume one line of the block."""
        if line.startswith("#"):
            self.flush()
            self.state = ""
        elif line.startswith("loop_"):
            self.flush()
            self.state = "head"
        elif line.startswith(("_chem_comp_atom.", "_chem_comp_bond.")) and self.state != "head":
            self.single_item(line)
        elif line.startswith("_chem_comp."):
            key, _, value = line.partition(" ")
            self.last = key.strip().removeprefix("_chem_comp.")
            self.meta[self.last] = value.strip()
            self.state, self.text = "meta", None
        elif self.state == "head" and line.startswith("_"):
            category, _, item = line.strip()[1:].partition(".")
            self.category = category
            self.items.append(item)
        elif self.state in ("head", "data"):
            self.loop_row(line)
        elif self.state == "meta":
            self.meta_text(line)

    def render(self) -> str:
        """Return the reduced block text."""
        self.flush()
        head = [self.out[0]] + [f"_chem_comp.{k} {self.meta.get(k, '?')}\n" for k in META]
        return "".join(head + self.out[1:]) + "#\n"


def main(source: str, snfg: str, out: str, release: str) -> None:
    """Write the reduced dictionary and report identifiers it lacks."""
    named = set(re.findall(r'"([A-Z0-9]{1,3})"', Path(snfg).read_text()))
    want = (set(CORE) | named) - DROP
    blocks: dict[str, Block] = {}
    current: Block | None = None
    with gzip.open(source, "rt") as handle:
        for line in handle:
            if line.startswith("data_"):
                name = line[5:].strip()
                current = Block(line) if name in want else None
                if current:
                    blocks[name] = current
            elif current:
                current.feed(line)
    body = "".join(blocks[name].render() for name in sorted(blocks))
    header = f"# molframe core components; wwPDB CCD release {release}\n"
    Path(out).write_text(header + body)
    absent = sorted(want - set(blocks))
    sys.stderr.write(f"{len(blocks)} blocks; not in dictionary: {absent}\n")


if __name__ == "__main__":
    main(*sys.argv[1:5])
