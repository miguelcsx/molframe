use pdbiox_core::structure::ChainSequenceExt;
use pdbiox_core::{InputBuffer, ReadOptions};

const ENTRY: &str = r"data_seq
_entity.id 1
_entity.type polymer
_struct_asym.id A
_struct_asym.entity_id 1
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 ALA
1 2 GLY
1 3 SER
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA ALA A X 1 0 0 0
2 C CA SER A X 3 1 0 0
_struct_ref.id R1
_struct_ref.entity_id 1
_struct_ref.db_name UNP
_struct_ref.db_code TEST_HUMAN
_struct_ref.pdbx_db_accession P00001
_struct_ref.pdbx_seq_one_letter_code 'A G S'
_struct_ref_seq.align_id 1
_struct_ref_seq.ref_id R1
_struct_ref_seq.pdbx_strand_id X
_struct_ref_seq.seq_align_beg 1
_struct_ref_seq.seq_align_end 3
_struct_ref_seq.db_align_beg 10
_struct_ref_seq.db_align_end 12
";

#[test]
fn chain_exposes_observed_canonical_reference_and_missing_positions() {
    let input = InputBuffer::from_bytes(ENTRY.as_bytes().to_vec());
    let structure = match crate::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let Some(chain) = structure.data().chains().next() else {
        panic!("chain absent")
    };
    assert_eq!(chain.observed_sequence().len(), 2);
    assert_eq!(chain.canonical_sequence().len(), 3);
    assert_eq!(
        chain.reference_sequences()[0].one_letter_code.as_deref(),
        Some("AGS")
    );
    assert_eq!(
        chain
            .sequence_mapping()
            .iter()
            .map(|mapping| mapping.reference_position)
            .collect::<Vec<_>>(),
        [Some(10), Some(12)]
    );
    assert_eq!(chain.missing_residues()[0].canonical_position, 2);
}

#[test]
fn canonical_writing_preserves_reference_sequences_and_alignments() {
    let input = InputBuffer::from_bytes(ENTRY.as_bytes().to_vec());
    let structure = match crate::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let canonical = crate::write_canonical(&structure);
    let input = InputBuffer::from_bytes(canonical.into_bytes());
    let round_tripped = match crate::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("canonical output failed: {findings:?}"),
    };
    let Some(chain) = round_tripped.data().chains().next() else {
        panic!("chain absent")
    };
    assert_eq!(
        chain.reference_sequences()[0].accession.as_deref(),
        Some("P00001")
    );
    assert_eq!(chain.sequence_mapping()[0].reference_position, Some(10));
}
