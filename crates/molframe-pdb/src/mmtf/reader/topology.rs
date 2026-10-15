//! Topology comparison between MMTF models.
//!
//! Complexity: O(atoms + residues + bonds) time per comparison; O(1) extra space.

use molframe_core::structure::StructureData;

/// Returns whether two structures have identical normalized topology.
///
/// The comparison covers declared hierarchy sizes, bonds, chain metadata,
/// residue identity fields, and atom identity fields. Coordinate-dependent
/// properties are intentionally excluded.
pub(super) fn same_topology(left: &StructureData, right: &StructureData) -> bool {
    if atom_count(left) != atom_count(right)
        || left.topology.residues.len() != right.topology.residues.len()
        || left.topology.chains.len() != right.topology.chains.len()
    {
        return false;
    }

    left.bonds.iter().eq(right.bonds.iter())
        && same_chain_topology(left, right)
        && same_residue_topology(left, right)
        && same_atom_topology(left, right)
}

/// Counts atoms without traversing the chunk columns.
pub(super) fn atom_count(data: &StructureData) -> u32 {
    data.chunks.last().map_or(0, |chunk| chunk.atoms().end)
}

/// Compares the topology-relevant metadata of every chain in two structures.
///
/// Chain labels, author labels, and polymer classifications must match in
/// iteration order.
pub(super) fn same_chain_topology(left: &StructureData, right: &StructureData) -> bool {
    left.chains().zip(right.chains()).all(|(left, right)| {
        left.label() == right.label()
            && left.auth_label() == right.auth_label()
            && left.polymer_kind() == right.polymer_kind()
    })
}

/// Compares the topology-relevant metadata of every residue in two structures.
///
/// Residue names, sequence identifiers, insertion codes, and hetero flags must
/// match in iteration order.
pub(super) fn same_residue_topology(left: &StructureData, right: &StructureData) -> bool {
    left.residues().zip(right.residues()).all(|(left, right)| {
        left.name() == right.name()
            && left.auth_seq_id() == right.auth_seq_id()
            && left.label_seq_id() == right.label_seq_id()
            && left.ins_code() == right.ins_code()
            && left.is_het() == right.is_het()
    })
}

/// Compares the topology-relevant metadata of every atom in two structures.
///
/// Atom names, elements, alternate labels, and formal charges must match in
/// iteration order.
pub(super) fn same_atom_topology(left: &StructureData, right: &StructureData) -> bool {
    left.atoms().zip(right.atoms()).all(|(left, right)| {
        left.name() == right.name()
            && left.element() == right.element()
            && left.alt_label() == right.alt_label()
            && left.formal_charge() == right.formal_charge()
    })
}
