//! Structure-comparison and sequence bindings.

mod compare;
mod compare_governed;
mod compare_metrics;
mod mapping;
mod sequence;
mod sequence_dispatch;
mod sequence_extended;
mod sequence_formats;

pub(crate) use compare::{
    PyCeAlignment, PyCeOptions, PyCeSignificanceProfile, PyDockQ, PyDockQOptions,
    PyEmptyLddtPolicy, PyEmptyQsPolicy, PyLddtOptions, PyQsOptions, analyse_dockq, analyse_gdt_ha,
    analyse_gdt_ts, analyse_lddt, analyse_qs_score, analyse_tm_score, ce_align, ce_alignments,
    dockq, gdt, gdt_ha, gdt_ts, lddt, qs_score, tm_score, weighted_rmsd,
};
pub(crate) use compare_governed::{
    analyse_ce_align, analyse_ce_alignments, analyse_gdt, analyse_map_chains, analyse_weighted_rmsd,
};
pub(crate) use compare_metrics::{
    PyCadContact, PyCadScore, PyContactArea, PyContactSimilarity, PyEquivalentAtomMapping,
    PyLigandRmsd, PyLocalCad, analyse_cad_contact_areas, analyse_cad_score,
    analyse_contact_map_similarity, analyse_equivalent_atom_mappings, analyse_ligand_symmetry_rmsd,
    cad_contact_areas, cad_score, contact_map_similarity, equivalent_atom_mappings,
    ligand_symmetry_rmsd,
};
pub(crate) use mapping::{
    PyChainAlternative, PyChainAssignment, PyChainMapping, PyChainSequence, PyResidueMatch,
    analyse_assign_chains, analyse_map_sequence_to_structure, assign_chains, chain_sequences,
    map_chains, map_sequence_to_structure,
};
pub(crate) use sequence::{PyAlignment, PyScoring, global, local, semi_global};
pub(crate) use sequence_dispatch::{
    PySequenceDocument, PySequenceDocumentKind, PySequenceFormat, read_sequence, write_sequence,
};
pub(crate) use sequence_extended::{
    PyFastaRecord, PyLadderDirection, PyMsaOptions, PySubstitutionMatrix, PyTree,
    a2m_match_columns, a3m_match_columns, align_global_banded, align_global_matrix,
    align_local_matrix, align_semi_global_matrix, blosum62, kmer_counts, minimizers,
    multiple_sequence_alignment, neighbor_joining, parse_a3m, parse_clustal, parse_fasta,
    parse_phylip, parse_stockholm, progressive_msa, seed_and_extend, upgma, write_a3m,
    write_clustal, write_fasta, write_phylip, write_stockholm,
};
pub(crate) use sequence_formats::{
    PyAlignmentMode, PyFastqRecord, PyMatrixProfile, PyRegionOptions, PySimilarKmer,
    PySimilarKmerOptions, align_region, load_matrix, parse_a2m, parse_fastq, similar_kmers,
    write_a2m, write_fastq,
};
