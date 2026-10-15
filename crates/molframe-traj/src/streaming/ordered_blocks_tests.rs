use super::*;

type Pipe = Ordered<u32, TrajectoryError>;

#[test]
fn blocks_are_delivered_in_order_whatever_the_completion_order() {
    let pipe = Pipe::new(3, 3, 1);
    assert_eq!(
        [pipe.claim(), pipe.claim(), pipe.claim(), pipe.claim()],
        [Some(0), Some(1), Some(2), None]
    );
    pipe.complete(2, Ok(22));
    pipe.complete(0, Ok(0));
    pipe.complete(1, Ok(11));
    assert_eq!(pipe.next(), Some(Ok(0)));
    assert_eq!(pipe.next(), Some(Ok(11)));
    assert_eq!(pipe.next(), Some(Ok(22)));
    assert_eq!(pipe.next(), None);
}

#[test]
fn claims_stay_inside_the_window_until_the_consumer_advances() {
    let pipe = Pipe::new(5, 2, 1);
    assert_eq!(pipe.claim(), Some(0));
    assert_eq!(pipe.claim(), Some(1));
    pipe.complete(0, Ok(7));
    assert_eq!(pipe.next(), Some(Ok(7)));
    assert_eq!(pipe.claim(), Some(2));
    pipe.stop();
    assert_eq!(pipe.claim(), None);
}

#[test]
fn a_block_error_ends_the_stream_after_earlier_blocks() {
    let pipe = Pipe::new(3, 3, 2);
    pipe.complete(1, Err(TrajectoryError::Cancelled));
    pipe.complete(0, Ok(5));
    assert_eq!(pipe.next(), Some(Ok(5)));
    assert_eq!(pipe.next(), Some(Err(TrajectoryError::Cancelled)));
    assert_eq!(pipe.next(), None);
}

#[test]
fn a_vanished_worker_is_reported_instead_of_waiting_forever() {
    let pipe = Pipe::new(2, 2, 1);
    pipe.retire();
    assert_eq!(pipe.next(), Some(Err(TrajectoryError::WorkerPanicked)));
    assert_eq!(pipe.next(), None);
}

#[test]
fn many_threads_deliver_every_block_exactly_once_in_order() {
    let pipe = Pipe::new(500, 6, 4);
    let delivered = std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                while let Some(id) = pipe.claim() {
                    pipe.complete(id, Ok(u32::try_from(id).unwrap_or(0)));
                }
                pipe.retire();
            });
        }
        let mut seen = Vec::new();
        while let Some(block) = pipe.next() {
            seen.push(block.expect("block"));
        }
        pipe.stop();
        seen
    });
    assert_eq!(delivered, (0..500).collect::<Vec<u32>>());
}
