use super::*;

#[test]
fn parentheses_operators_and_quoted_values_are_distinct_tokens() {
    let tokens = match lex("(name 'C A') and occupancy >= 0.5") {
        Ok(tokens) => tokens,
        Err(error) => panic!("lex failed: {error}"),
    };
    assert_eq!(tokens.len(), 8);
    assert!(matches!(tokens[0].kind, TokenKind::LeftParen));
    assert_eq!(tokens[2].kind, TokenKind::Value("C A".into()));
    assert_eq!(tokens[6].kind, TokenKind::Operator(">=".into()));
}

#[test]
fn an_unclosed_quote_is_a_selection_diagnostic() {
    let result = lex("name 'CA");
    assert_eq!(
        result.err().map(|finding| finding.code()),
        Some(Code::E4001)
    );
}

#[test]
fn every_token_carries_its_byte_range_line_and_column() {
    let tokens = match lex("chain A\n  and (name 'C A')") {
        Ok(tokens) => tokens,
        Err(error) => panic!("lex failed: {error}"),
    };
    let spans = tokens
        .iter()
        .map(|token| {
            (
                token.span.start.byte_offset,
                token.span.end,
                token.span.start.line,
                token.span.start.column,
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        spans,
        vec![
            (0, 5, 1, 1),
            (6, 7, 1, 7),
            (10, 13, 2, 3),
            (14, 15, 2, 7),
            (15, 19, 2, 8),
            (20, 25, 2, 13),
            (25, 26, 2, 18),
        ]
    );
}

#[test]
fn an_unclosed_quote_spans_from_the_quote_to_the_end() {
    let Err(finding) = lex("name 'CA") else {
        panic!("an unclosed quote is an error");
    };
    let Some(span) = finding.span() else {
        panic!("the finding is located");
    };
    assert_eq!((span.start.byte_offset, span.end), (5, 8));
}
