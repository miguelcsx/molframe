use super::*;
use crate::msa::profile::GAP;

fn scoring() -> Scoring {
    Scoring {
        match_score: 2,
        mismatch_score: -1,
        gap_open: -3,
        gap_extend: -1,
    }
}

#[test]
fn progressive_alignment_preserves_input_order_and_ungapped_sequences() {
    let sequences: [&[u8]; 4] = [b"ACGT", b"ACGGT", b"AGT", b"TTTT"];
    let result = progressive_msa(
        &sequences,
        MsaOptions::progressive(scoring()).with_refinement_passes(2),
    );
    let Ok(rows) = result else {
        panic!("valid progressive MSA");
    };
    assert!(rows.windows(2).all(|pair| pair[0].len() == pair[1].len()));
    for (row, original) in rows.iter().zip(sequences) {
        let ungapped = row
            .iter()
            .copied()
            .filter(|symbol| *symbol != GAP)
            .collect::<Vec<_>>();
        assert_eq!(ungapped, original);
    }
}

#[test]
fn guide_tree_merges_the_similar_sequences_before_the_outgroup() {
    let sequences: [&[u8]; 3] = [b"AAAA", b"AAAT", b"TTTT"];
    let Ok(distances) = pairwise_distances(&sequences, scoring()) else {
        panic!("small fixture dimensions fit");
    };
    let members = vec![Some(vec![0]), Some(vec![1]), Some(vec![2])];
    assert_eq!(closest_pair(&members, &distances), Ok(Some((0, 1))));
}

#[test]
fn refinement_never_reduces_the_sum_of_pairs_objective() {
    let sequences: [&[u8]; 4] = [b"ABCD", b"ABXCD", b"ACD", b"ABYD"];
    let base = progressive_msa(&sequences, MsaOptions::progressive(scoring()));
    let refined = progressive_msa(
        &sequences,
        MsaOptions::progressive(scoring()).with_refinement_passes(3),
    );
    let (Ok(base), Ok(refined)) = (base, refined) else {
        panic!("valid MSA");
    };
    let base_profile = Profile {
        rows: base.into_iter().enumerate().collect(),
    };
    let refined_profile = Profile {
        rows: refined.into_iter().enumerate().collect(),
    };
    let (Ok(refined_score), Ok(base_score)) = (
        refined_profile.score(scoring()),
        base_profile.score(scoring()),
    ) else {
        panic!("fixture score is representable")
    };
    assert!(refined_score >= base_score);
}

#[test]
fn a_gap_reward_is_rejected_explicitly() {
    let invalid = Scoring {
        gap_open: 1,
        ..scoring()
    };
    assert_eq!(
        progressive_msa(&[b"A", b"A"], MsaOptions::progressive(invalid)),
        Err(MsaError::InvalidGapScore)
    );
}

#[test]
fn empty_and_single_inputs_have_total_results() {
    assert_eq!(
        progressive_msa(&[], MsaOptions::progressive(scoring())),
        Ok(Vec::new())
    );
    assert_eq!(
        progressive_msa(&[b"ABC"], MsaOptions::progressive(scoring())),
        Ok(vec![b"ABC".to_vec()])
    );
}

#[test]
fn a_deletion_and_an_insertion_open_shared_profile_columns() {
    let deletion = msa(&[b"ACGT", b"AGT", b"ACGT"], scoring());
    let Ok(deletion) = deletion else {
        panic!("valid deletion MSA");
    };
    assert_eq!(deletion[1], b"A-GT");
    assert!(deletion.iter().all(|row| row.len() == 4));

    let insertion = msa(&[b"ACGT", b"ACXGT"], scoring());
    let Ok(insertion) = insertion else {
        panic!("valid insertion MSA");
    };
    assert_eq!(insertion[0], b"AC-GT");
    assert_eq!(insertion[1], b"ACXGT");
}
