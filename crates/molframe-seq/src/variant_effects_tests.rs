use super::*;
use crate::variants::parse_variants;

const TABLE: &str =
    "mutant,DMS_score,pred\nA1G,-0.5,0.1\nA1V,-0.9,0.2\nG2A,0.3,\nA1G:G2A,-1.2,0.5\nfoo,1,1\n";

#[test]
fn positions_are_summed_from_single_substitutions_only() {
    let t = parse_variants(TABLE, "t.csv").expect("a result");
    let e = per_position(&t, 0);
    assert_eq!(e.len(), 2);
    assert_eq!((e[0].position, e[0].count), (1, 2));
    assert!((e[0].mean + 0.7).abs() < 1e-9);
    assert_eq!((e[0].min, e[0].max), (-0.9, -0.5));
    assert_eq!(
        e[1].count, 1,
        "the double mutant is attributable to neither position"
    );
    let pred = per_position(&t, 1);
    assert_eq!(
        pred.iter().find(|x| x.position == 2).map(|x| x.count),
        None,
        "G2A has no prediction"
    );
}

#[test]
fn the_numbering_is_checked_against_the_structure_and_an_offset_is_found_when_it_differs() {
    let t = parse_variants(TABLE, "t.csv").expect("a result");
    let e = per_position(&t, 0);
    // the structure numbers the same two residues 11 and 12
    let residues = [(11, 'A'), (12, 'G'), (13, 'L')];
    let none = map_positions(&e, &residues, 0);
    assert_eq!((none.agree, none.absent), (0, 2));
    let best = best_offset(&e, &residues, 50).expect("a result");
    assert_eq!(best.offset, 10);
    assert_eq!(best.agreement(), Some(1.0));
    // a wrong letter is counted as a disagreement, not forced to fit
    let wrong = map_positions(&e, &[(11, 'A'), (12, 'W')], 10);
    assert_eq!((wrong.agree, wrong.disagree), (1, 1));
    assert_eq!(wrong.agreement(), Some(0.5));
    assert!(best_offset(&e, &[(11, 'W'), (12, 'W')], 50).is_none());
}

#[test]
fn an_offset_is_found_across_the_whole_range_in_either_direction() {
    let t = parse_variants("mutant,s\nA300G,1\nG301A,1\n", "t.csv").expect("a result");
    let e = per_position(&t, 0);
    // the structure numbers the same residues 74 and 75: offset -226
    let residues = [(74, 'A'), (75, 'G')];
    assert_eq!(
        best_offset(&e, &residues, 226).expect("a result").offset,
        -226
    );
    assert!(
        best_offset(&e, &residues, 225).is_none(),
        "range means plus or minus range, no more"
    );
    let up = [(500, 'A'), (501, 'G')];
    assert_eq!(best_offset(&e, &up, 200).expect("a result").offset, 200);
}

#[test]
fn a_position_carries_its_median_and_the_share_past_the_stated_cutoff() {
    let t = parse_variants(
        "mutant,s\nA1G,-1\nA1V,0.5\nA1L,0.7\nA1D,-3\nA1K,0.2\n",
        "t.csv",
    )
    .expect("a result");
    let e = per_position(&t, 0);
    assert!((e[0].median - 0.2).abs() < 1e-9);
    assert!((e[0].fraction_deleterious - 0.4).abs() < 1e-9);
    let high = per_position_at(
        &t,
        0,
        Cutoff {
            threshold: 0.5,
            low_is_deleterious: false,
        },
    );
    assert!((high[0].fraction_deleterious - 0.4).abs() < 1e-9);
    assert_eq!(Cutoff::default().text(), "deleterious = value <= 0");
    let even = parse_variants("mutant,s\nA1G,1\nA1V,3\n", "t.csv").expect("a result");
    assert!((per_position(&even, 0)[0].median - 2.0).abs() < 1e-9);
}
