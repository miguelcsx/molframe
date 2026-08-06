use super::*;

#[test]
fn interning_the_same_string_twice_yields_the_same_identifier() -> Result<(), DictionaryFull> {
    let mut interner = Interner::new();
    let first = interner.intern("MY_LIGAND")?;
    let second = interner.intern("MY_LIGAND")?;
    assert_eq!(first, second);
    assert_eq!(interner.len(), 1);
    Ok(())
}

#[test]
fn canonical_strings_never_enter_the_local_dictionary() -> Result<(), DictionaryFull> {
    let mut interner = Interner::new();
    for text in ["ALA", "CA", "HOH", "ZN"] {
        let id = interner.intern(text)?;
        assert!(id.is_canonical(), "{text} should be canonical");
    }
    assert!(interner.is_empty());
    Ok(())
}

#[test]
fn two_interners_agree_on_canonical_identifiers_and_not_on_local_ones() -> Result<(), DictionaryFull>
{
    let (mut left, mut right) = (Interner::new(), Interner::new());
    assert_eq!(left.intern("GLY")?, right.intern("GLY")?);

    let left_ligand = left.intern("LIG")?;
    right.intern("OTHER")?;
    let right_ligand = right.intern("LIG")?;
    assert_ne!(left_ligand, right_ligand);
    Ok(())
}

#[test]
fn local_identifiers_are_issued_in_first_seen_order() -> Result<(), DictionaryFull> {
    let mut interner = Interner::new();
    let first = interner.intern("AAA")?;
    let second = interner.intern("BBB")?;
    interner.intern("AAA")?;
    let third = interner.intern("CCC")?;
    assert!(first < second && second < third);
    Ok(())
}

#[test]
fn every_local_identifier_resolves_to_the_string_it_was_made_from() -> Result<(), DictionaryFull> {
    let mut interner = Interner::new();
    let names = ["O5'", "MY_LAB_ATOM", "1HB", "", "a very long ligand code"];
    let ids: Vec<_> = names
        .iter()
        .map(|name| interner.intern(name))
        .collect::<Result<_, _>>()?;
    for (id, name) in ids.iter().zip(names) {
        assert_eq!(interner.resolve(*id), Some(name));
    }
    Ok(())
}

#[test]
fn a_string_not_yet_seen_resolves_to_nothing_without_growing_the_dictionary() {
    let interner = Interner::new();
    assert_eq!(interner.get("NEVER_SEEN"), None);
    assert!(interner.is_empty());
    assert!(interner.get("ALA").is_some());
}

#[test]
fn the_dictionary_refuses_to_grow_past_its_limit_rather_than_exhausting_memory()
-> Result<(), DictionaryFull> {
    let mut interner = Interner::new().with_limit(2);
    interner.intern("one")?;
    interner.intern("two")?;
    assert_eq!(interner.intern("three"), Err(DictionaryFull));
    // A string already present is still resolvable once the limit is reached.
    assert!(interner.intern("one").is_ok());
    Ok(())
}

#[test]
fn the_arena_holds_local_strings_once_rather_than_one_allocation_each() -> Result<(), DictionaryFull>
{
    let mut interner = Interner::new();
    interner.intern("abc")?;
    interner.intern("de")?;
    interner.intern("abc")?;
    assert_eq!(interner.arena_len(), 5);
    Ok(())
}
