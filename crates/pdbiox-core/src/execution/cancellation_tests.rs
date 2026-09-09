use super::*;

#[test]
fn cancellation_is_visible_through_every_clone() {
    let token = CancellationToken::new();
    let clone = token.clone();
    assert_eq!(clone.check(), Ok(()));
    token.cancel();
    assert_eq!(clone.check(), Err(Cancelled));
}
