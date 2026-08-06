use super::*;
use crate::symbol::{DictionaryFull, Interner};

#[test]
fn the_blank_altloc_is_distinct_from_every_label() -> Result<(), DictionaryFull> {
    let mut interner = Interner::new();
    let a = AltId::labelled(interner.intern("A")?);
    assert!(AltId::BLANK.is_blank());
    assert!(!a.is_blank());
    assert_eq!(AltId::BLANK.symbol(), None);
    assert_eq!(a.symbol(), Some(interner.intern("A")?));
    Ok(())
}

#[test]
fn an_altloc_label_longer_than_one_character_is_representable() -> Result<(), DictionaryFull> {
    let mut interner = Interner::new();
    let long = AltId::labelled(interner.intern("ALT2")?);
    let resolved = long.symbol().and_then(|symbol| interner.resolve(symbol));
    assert_eq!(resolved, Some("ALT2"));
    Ok(())
}
