use crate::{HeteroPolicy, MetalFallback};

#[test]
fn named_and_listed_hetero_classes_round_trip() {
    for text in ["none", "solute", "all", "waters", "waters,ions", "ligands"] {
        let Ok(policy) = text.parse::<HeteroPolicy>() else {
            panic!("`{text}` is a valid hetero policy");
        };
        assert_eq!(policy.to_string(), text);
    }
}

#[test]
fn hetero_words_combine_as_a_union() {
    let Ok(policy) = "solute, waters".parse::<HeteroPolicy>() else {
        panic!("valid");
    };
    assert_eq!(policy, HeteroPolicy::ALL);
}

#[test]
fn an_unknown_word_is_refused_with_the_choices() {
    let Err(error) = "lipids".parse::<HeteroPolicy>() else {
        panic!("`lipids` is not a class");
    };
    assert!(error.to_string().contains("waters, ligands, ions"));
    assert!("copper".parse::<MetalFallback>().is_err());
}

#[test]
fn metal_fallback_names_round_trip() {
    for fallback in [
        MetalFallback::None,
        MetalFallback::Vdw,
        MetalFallback::Ionic,
    ] {
        assert_eq!(fallback.to_string().parse(), Ok(fallback));
    }
}
