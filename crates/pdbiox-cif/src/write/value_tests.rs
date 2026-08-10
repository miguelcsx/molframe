use super::*;
use crate::{CifValue, Document, parse};
use pdbiox_core::io::InputBuffer;

#[test]
fn every_quote_shape_reparses_to_the_original_text() {
    for text in [
        "plain",
        "two words",
        "don't",
        "say \"don't\"",
        "line\nbreak",
    ] {
        let rendered = render_value(&CifValue::Text(text.to_owned().into()));
        let source = format!("data_x\n_test.value {rendered}\n");
        let input = InputBuffer::from_bytes(source.into_bytes());
        let Ok((document, _)) = parse(&input) else {
            panic!("rendered value did not parse: {rendered}");
        };
        assert_eq!(value(&document), Some(text));
    }
}

#[test]
fn reserved_cif_words_are_not_rendered_bare() {
    for text in ["loop_", "STOP_", "data_x", "save_frame", ".", "?"] {
        assert_ne!(quote_text(text), text);
    }
}

fn value(document: &Document) -> Option<&str> {
    document.first_block()?.category("test")?.text("value", 0)
}
