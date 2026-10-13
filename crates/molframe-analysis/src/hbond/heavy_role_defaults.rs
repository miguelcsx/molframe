//! The built-in donor and acceptor typing for heavy-atom hydrogen bonds, as data.
//!
//! The rows follow the conventions of `McDonald` and Thornton's `HBPLUS` and of
//! Baker and Hubbard. Nothing about residues or atoms is decided anywhere else:
//! the typing code reads a [`HeavyRoleTable`](super::HeavyRoleTable) and this is
//! the table it reads unless the caller supplies another, in the same text
//! format, to [`HeavyRoleTable::parse`](super::HeavyRoleTable::parse).

use molframe_chem::TableIdentity;

/// Identity of the built-in table.
pub const HEAVY_ROLES_DEFAULT: TableIdentity = TableIdentity::builtin(
    "hbplus-baker-hubbard-heavy-roles",
    "1",
    "McDonald IK, Thornton JM (1994) J Mol Biol 238:777 (HBPLUS) and Baker EN, Hubbard RE (1984) Prog Biophys Mol Biol 44:97, as encoded in this crate",
    "MIT (this repository)",
);

/// The built-in rows, in the text format a caller's own table uses.
///
/// Protonation variants are matched on their own name before their parent's; for one
/// residue and atom the first row whose condition holds applies. A `weak` row is
/// skipped unless sulfur is admitted, and a `chain-start` row only applies to the
/// first residue of a chain.
pub const HEAVY_ROLES_DEFAULT_TEXT: &str = "\
amino: ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
nucleic: DA A DC C DG G DT DU U DI I
water: HOH WAT H2O DOD

# residues      atoms        role      [condition]
# protonation variants
HID,HSD         ND1          donor
HID,HSD         NE2          acceptor
HIE,HSE         NE2          donor
HIE,HSE         ND1          acceptor
HIP,HSP         ND1,NE2      donor
ASH             OD2          both
GLH             OE2          both
LYN             NZ           both
ASH             OD1          acceptor
GLH             OE1          acceptor
CYX,CYM         SG           acceptor  weak

# water: an oxygen does both
@water          O*           both

# protein backbone and side chains
PRO             N            donor     chain-start
PRO             N            none
@amino          N            donor
ARG             NE,NH1,NH2   donor
ASN             ND2          donor
GLN             NE2          donor
LYS             NZ           donor
TRP             NE1          donor
@amino          O,OXT        acceptor
ASN             OD1          acceptor
GLN             OE1          acceptor
ASP             OD1,OD2      acceptor
GLU             OE1,OE2      acceptor
HIS             ND1,NE2      both
SER             OG           both
THR             OG1          both
TYR             OH           both
CYS             SG           both      weak
MET             SD           acceptor  weak

# nucleic acids: phosphate and ester oxygens, then the bases
@nucleic        OP1,OP2,OP3,O1P,O2P,O3P,O5',O3',O4'  acceptor
@nucleic        O2'          both
DA,A            N6           donor
DC,C            N4           donor
DG,G            N2           donor
DG,G,DI,I       N1           donor
DT,DU,U         N3           donor
DA,A            N1,N3,N7     acceptor
DC,C            N3,O2        acceptor
DG,G,DI,I       N3,N7,O6     acceptor
DT,DU,U         O2,O4        acceptor
";
