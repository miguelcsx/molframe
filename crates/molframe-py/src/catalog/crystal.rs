//! Catalogued operations of the `crystal` namespace.

use super::{Capability, Domain};

const D: Domain = Domain("crystal");

pub(super) const ROWS: &[Capability] = &[
    D.row(
        "assembly_covalent_links",
        "structure,assembly,model",
        "AssemblyBond[]",
        "cpu,spatial,memory",
    ),
    D.row("assemblies", "structure", "str[]", "cpu"),
    D.row(
        "assembly",
        "structure,id",
        "AssemblyInstance[]",
        "cpu,memory",
    ),
    D.row(
        "normalizers",
        "cell,space_group,hkl,amplitudes,bins",
        "Float64Array",
        "cpu,memory",
    ),
    D.row(
        "reduce_cell",
        "lengths,angles,epsilon,iteration_limit",
        "ReducedCell",
        "cpu",
    ),
    D.row(
        "structure_factors",
        "structure,hkl",
        "Complex128Array",
        "cpu,memory",
    ),
    D.row("read_mrc", "path or bytes", "DensityMap", "cpu,memory"),
    D.row("read_mtz", "path or bytes", "ReflectionTable", "cpu,memory"),
    D.row("space_group", "structure", "SpaceGroup | None", "cpu"),
    D.row(
        "assembly_chunks",
        "structure,assembly",
        "AssemblyChunks",
        "cpu,memory",
    ),
    D.row(
        "assembly_structure",
        "structure,assembly",
        "Structure",
        "cpu,memory",
    ),
    D.row("ncs_chain_pairs", "structure", "NcsChainPair[]", "cpu"),
    D.row("ncs_operators", "structure", "NcsOperator[]", "cpu"),
    D.row("ncs_structure", "structure", "Structure", "cpu,memory"),
];
