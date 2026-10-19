//! Catalogued operations of the `chemistry` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("chemistry");

pub(super) const ROWS: &[Capability] = &[
    D.row(
        "apply_polymer_role_profile",
        "structure,components,rules,profile_id,version",
        "PolymerRoleReport",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "partial_charges",
        "structure,components?",
        "PartialCharges",
        "cpu,memory",
    ),
    D.row("element", "symbol", "ElementProperties", "cpu"),
    D.row("vdw_radius", "symbol", "float", "cpu"),
    D.row("vdw_radii", "structure", "array", "cpu"),
    D.row("covalent_radius", "symbol", "float", "cpu"),
    D.row("covalent_radii", "structure", "array", "cpu"),
    D.row(
        "carbohydrates",
        "structure, component dictionary",
        "CarbohydrateReport",
        "cpu,spatial,memory",
    ),
    D.row(
        "snfg_symbol",
        "component identifier",
        "SnfgSymbol | None",
        "cpu",
    ),
    D.row(
        "annotate",
        "structure,components,version,polymer_links?",
        "Structure",
        "cpu,memory",
    )
    .policy(),
    D.row("polymer_atom_roles", "none", "dict[str,int]", "cpu"),
    D.row("read_sdf", "path or text", "list[Molecule]", "cpu"),
    D.row("read_mol", "path or text", "Molecule", "cpu"),
    D.row("write_sdf", "molecules", "text", "cpu"),
    D.row("molecule", "structure", "Molecule", "cpu"),
    D.row("smarts", "structure,pattern", "list[tuple]", "cpu"),
];
