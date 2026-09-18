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
