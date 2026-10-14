use super::*;

const ROW: Row = Row::new("lab table", 7);

#[test]
fn words_are_counted_exactly() {
    assert_eq!(ROW.words::<2>("C 1.7", "`a b`").ok(), Some(["C", "1.7"]));
    for text in ["C", "C 1 2", ""] {
        let Err(error) = ROW.words::<2>(text, "`a b`") else {
            panic!("{text:?} accepted");
        };
        assert_eq!(error.line(), 7);
    }
}

#[test]
fn numbers_are_finite_and_bounded() {
    assert_eq!(ROW.finite::<f32>("1.5").ok(), Some(1.5));
    assert!(ROW.finite::<f64>("nan").is_err());
    assert!(ROW.finite::<f32>("x").is_err());
    assert_eq!(ROW.bounded::<f64>("2", 3.0, "area").ok(), Some(2.0));
    for word in ["0", "-1", "4", "inf"] {
        assert!(ROW.bounded::<f64>(word, 3.0, "area").is_err(), "{word}");
    }
}

#[test]
fn elements_resolve_and_keys_are_claimed_once() {
    assert!(ROW.element("Fe").is_ok());
    assert!(ROW.element("Qq").is_err());
    let mut slots = [None, None];
    assert!(ROW.claim(slots.get_mut(0), "a", 1).is_ok());
    assert!(ROW.claim(slots.get_mut(0), "a", 2).is_err());
    assert!(ROW.claim(slots.get_mut(5), "z", 2).is_err());
    assert_eq!(slots, [Some(1), None]);
}
