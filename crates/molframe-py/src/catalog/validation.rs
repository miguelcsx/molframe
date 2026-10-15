//! Catalogued operations of the `validation` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("validation");

pub(super) const ROWS: &[Capability] = &[
    D.row("clashes", "structure", "ClashTable", "cpu"),
    D.row(
        "bond_length_deviations",
        "structure,tolerance",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "cis_peptides",
        "structure,threshold_degrees",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "quality_flags",
        "structure",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row("valence", "structure", "Analysis<Table>", "cpu,memory")
        .policy(),
    D.row(
        "planarity",
        "structure,max_deviation",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row("completeness", "structure", "Analysis<Table>", "cpu,memory")
        .policy(),
    D.row(
        "altloc_occupancy_sums",
        "structure,expected_sum,tolerance",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "ligand_geometry",
        "structure,tolerance",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "b_factor_distribution",
        "structure,outlier_standard_deviations,selection",
        "Analysis<Table>",
        "cpu,memory",
    )
    .policy(),
    D.row(
        "clashscore",
        "structure",
        "ValidationReport",
        "cpu,spatial,memory",
    ),
    D.row(
        "geometry_z",
        "structure",
        "ValidationReport",
        "cpu,spatial,memory",
    ),
    D.row(
        "nucleic_conformation",
        "structure",
        "ValidationReport",
        "cpu,spatial,memory",
    ),
    D.row(
        "protein_conformation",
        "structure",
        "ValidationReport",
        "cpu,spatial,memory",
    ),
];
