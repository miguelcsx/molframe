use super::*;

fn tokens(text: &str) -> Vec<Token<'_>> {
    let mut lexer = match Lexer::new(text.as_bytes()) {
        Ok(lexer) => lexer,
        Err(error) => panic!("lexing failed: {error:?}"),
    };
    let mut found = Vec::new();
    while let Ok(Some(spanned)) = lexer.next_token() {
        found.push(spanned.token);
    }
    found
}

#[test]
fn a_block_header_a_tag_and_a_value_are_three_tokens() {
    assert_eq!(
        tokens("data_test\n_entry.id 1ABC\n"),
        [
            Token::Block("test"),
            Token::Tag("_entry.id"),
            Token::Value("1ABC", Quoting::Bare)
        ]
    );
}

#[test]
fn comments_and_blank_lines_are_skipped_between_tokens() {
    assert_eq!(
        tokens("# a comment\n\ndata_x\n  # another\n_a.b 1\n"),
        [
            Token::Block("x"),
            Token::Tag("_a.b"),
            Token::Value("1", Quoting::Bare)
        ]
    );
}

#[test]
fn quoting_travels_with_the_value_so_a_sentinel_can_be_told_from_a_string() {
    assert_eq!(
        tokens("data_x\n_a.b . _a.c '.' _a.d \"?\"\n"),
        [
            Token::Block("x"),
            Token::Tag("_a.b"),
            Token::Value(".", Quoting::Bare),
            Token::Tag("_a.c"),
            Token::Value(".", Quoting::Single),
            Token::Tag("_a.d"),
            Token::Value("?", Quoting::Double),
        ]
    );
}

#[test]
fn a_quote_inside_a_word_does_not_close_the_value() {
    assert_eq!(
        tokens("data_x\n_a.b 'it's here'\n"),
        [
            Token::Block("x"),
            Token::Tag("_a.b"),
            Token::Value("it's here", Quoting::Single)
        ]
    );
}

#[test]
fn a_semicolon_block_carries_its_lines_and_ends_at_the_next_semicolon() {
    let found = tokens("data_x\n_a.b\n;first line\nsecond line\n;\n");
    assert_eq!(
        found.get(2),
        Some(&Token::Value("first line\nsecond line", Quoting::Text))
    );
}

#[test]
fn a_semicolon_that_is_not_in_the_first_column_is_ordinary_text() {
    assert_eq!(
        tokens("data_x\n_a.b a;b\n"),
        [
            Token::Block("x"),
            Token::Tag("_a.b"),
            Token::Value("a;b", Quoting::Bare)
        ]
    );
}

#[test]
fn loop_and_save_frames_are_recognised_whatever_their_case() {
    assert_eq!(
        tokens("DATA_x\nLOOP_\n_a.b\n1\nsave_frame\nsave_\n"),
        [
            Token::Block("x"),
            Token::Loop,
            Token::Tag("_a.b"),
            Token::Value("1", Quoting::Bare),
            Token::FrameStart("frame"),
            Token::FrameEnd,
        ]
    );
}

#[test]
fn a_quote_that_never_closes_is_reported_with_where_it_opened() {
    let mut lexer = match Lexer::new(b"data_x\n_a.b 'unclosed\n") {
        Ok(lexer) => lexer,
        Err(error) => panic!("{error:?}"),
    };
    let _ = lexer.next_token();
    let _ = lexer.next_token();
    assert!(matches!(
        lexer.next_token(),
        Err(LexError::UnterminatedQuote(_))
    ));
}

#[test]
fn a_text_block_that_never_closes_is_reported() {
    let mut lexer = match Lexer::new(b"data_x\n_a.b\n;forever\n") {
        Ok(lexer) => lexer,
        Err(error) => panic!("{error:?}"),
    };
    let _ = lexer.next_token();
    let _ = lexer.next_token();
    assert!(matches!(
        lexer.next_token(),
        Err(LexError::UnterminatedText(_))
    ));
}

#[test]
fn every_token_carries_the_line_it_was_found_on() {
    let mut lexer = match Lexer::new(b"data_x\n\n_a.b 1\n") {
        Ok(lexer) => lexer,
        Err(error) => panic!("{error:?}"),
    };
    let first = lexer
        .next_token()
        .ok()
        .flatten()
        .map(|spanned| spanned.span.start.line);
    let second = lexer
        .next_token()
        .ok()
        .flatten()
        .map(|spanned| spanned.span.start.line);
    assert_eq!((first, second), (Some(1), Some(3)));
}

#[test]
fn positions_are_exact_for_lf_crlf_and_multibyte_bare_words() {
    let mut lexer = match Lexer::new("data_α\r\n_entry.id β\r\nloop_\n_a.x 1\n".as_bytes()) {
        Ok(lexer) => lexer,
        Err(error) => panic!("{error:?}"),
    };
    let expected = [
        Spanned {
            token: Token::Block("α"),
            span: ByteSpan::new(Position::new(0, 1, 1), 7),
        },
        Spanned {
            token: Token::Tag("_entry.id"),
            span: ByteSpan::new(Position::new(9, 2, 1), 18),
        },
        Spanned {
            token: Token::Value("β", Quoting::Bare),
            span: ByteSpan::new(Position::new(19, 2, 11), 21),
        },
        Spanned {
            token: Token::Loop,
            span: ByteSpan::new(Position::new(23, 3, 1), 28),
        },
        Spanned {
            token: Token::Tag("_a.x"),
            span: ByteSpan::new(Position::new(29, 4, 1), 33),
        },
        Spanned {
            token: Token::Value("1", Quoting::Bare),
            span: ByteSpan::new(Position::new(34, 4, 6), 35),
        },
    ];
    for expected in expected {
        assert_eq!(lexer.next_token(), Ok(Some(expected)));
    }
    assert_eq!(lexer.next_token(), Ok(None));
    assert_eq!(lexer.position(), Position::new(36, 5, 1));
}

#[test]
fn input_that_is_not_text_is_refused_rather_than_guessed_at() {
    assert!(matches!(Lexer::new(&[0xff, 0xfe]), Err(LexError::NotText)));
}
