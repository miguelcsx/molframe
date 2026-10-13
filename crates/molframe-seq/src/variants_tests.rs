use super::*;

#[test]
fn substitutions_are_read_in_the_forms_tables_use() {
    let s = |t: &str| parse_substitution(t);
    assert_eq!(
        s("A123G"),
        Some(Substitution {
            wild_type: 'A',
            position: 123,
            mutant: 'G'
        })
    );
    assert_eq!(
        s("p.Ala123Gly"),
        Some(Substitution {
            wild_type: 'A',
            position: 123,
            mutant: 'G'
        })
    );
    assert_eq!(s("K5*").expect("a stop").mutant, '*');
    assert_eq!(s("Lys5Ter").expect("a stop").mutant, '*');
    assert_eq!(s("a7v").expect("lower case").wild_type, 'A');
    assert!(
        s("hello").is_none() && s("A12").is_none() && s("123").is_none() && s("AB12G").is_none()
    );
}

const TABLE: &str =
    "mutant,DMS_score,pred\nA1G,-0.5,0.1\nA1V,-0.9,0.2\nG2A,0.3,\nA1G:G2A,-1.2,0.5\nfoo,1,1\n";

#[test]
fn a_table_gives_columns_variants_and_says_what_it_left_out() {
    let t = parse_variants(TABLE, "t.csv").expect("a readable table");
    assert_eq!(t.columns, ["DMS_score", "pred"]);
    assert_eq!(t.variants.len(), 4);
    assert_eq!(t.variants[3].changes.len(), 2);
    assert_eq!(
        t.variants[2].values,
        [Some(0.3), None],
        "an empty cell is absent, not zero"
    );
    assert!(t.notes.iter().any(|n| n.contains("foo")));
}

#[test]
fn a_table_that_cannot_be_read_says_why() {
    assert!(
        parse_variants("", "x")
            .expect_err("unreadable")
            .to_string()
            .contains("empty")
    );
    assert!(
        parse_variants("name,score\na,1\n", "x")
            .expect_err("unreadable")
            .to_string()
            .contains("substitutions")
    );
    assert!(
        parse_variants("mutant,note\nA1G,hello\n", "x")
            .expect_err("unreadable")
            .to_string()
            .contains("numeric")
    );
    assert!(
        parse_variants("mutant,score\nnot,1\n", "x")
            .expect_err("unreadable")
            .to_string()
            .contains("A123G")
    );
}

#[test]
fn quoted_cells_may_hold_the_delimiter_and_cells_are_trimmed() {
    assert_eq!(
        split_delimited_row(r#" a ,"b,c","say ""hi""", d"#, ','),
        ["a", "b,c", r#"say "hi""#, "d"]
    );
}
