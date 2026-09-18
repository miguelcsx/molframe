use super::*;
use crate::{Sample, structure};

fn source(copies: u64) -> SyntheticCifSource {
    let tile = Tile::from_structure(&structure(Sample::Tiny), Seed::new(1));
    SyntheticCifSource::new(tile, copies)
}

fn drain(mut stream: SyntheticCifSource) -> String {
    let mut text = String::new();
    if let Err(error) = stream.read_to_string(&mut text) {
        panic!("reading the synthetic stream failed: {error}")
    }
    text
}

#[test]
fn the_stream_opens_with_a_block_header_and_a_loop() {
    let text = drain(source(1));
    assert!(text.starts_with("data_synthetic\nloop_\n"), "bad header");
    assert!(text.contains("_atom_site.Cartn_x\n"));
}

#[test]
fn one_coordinate_row_is_emitted_per_atom_of_every_copy() {
    let tile = Tile::from_structure(&structure(Sample::Tiny), Seed::new(1));
    let expected = tile.atoms().len() * 3;
    let text = drain(SyntheticCifSource::new(tile, 3));
    let rows = text
        .lines()
        .filter(|line| line.starts_with("ATOM "))
        .count();
    assert_eq!(rows, expected);
}

#[test]
fn every_row_carries_one_field_per_declared_item() {
    let text = drain(source(1));
    let items = text
        .lines()
        .filter(|line| line.starts_with("_atom_site."))
        .count();
    for line in text.lines().filter(|line| line.starts_with("ATOM ")) {
        assert_eq!(
            line.split_whitespace().count(),
            items,
            "field count differs from the header: {line}"
        );
    }
}

#[test]
fn serial_numbers_ascend_without_a_gap() {
    let text = drain(source(2));
    for (expected, line) in (1_u64..).zip(text.lines().filter(|line| line.starts_with("ATOM "))) {
        let Some(field) = line.split_whitespace().nth(1) else {
            panic!("row has no serial: {line}")
        };
        let Ok(serial) = field.parse::<u64>() else {
            panic!("serial is not a number: {field}")
        };
        assert_eq!(serial, expected);
    }
}

#[test]
fn the_reported_atom_count_matches_what_is_emitted() {
    let stream = source(4);
    let declared = stream.atom_count();
    let rows = drain(stream)
        .lines()
        .filter(|line| line.starts_with("ATOM "))
        .count();
    let Ok(rows) = u64::try_from(rows) else {
        panic!("row count does not fit a u64")
    };
    assert_eq!(declared, rows);
}

#[test]
fn the_same_seed_produces_the_same_bytes() {
    assert_eq!(drain(source(2)), drain(source(2)));
}

#[test]
fn a_different_seed_produces_different_bytes() {
    let other = SyntheticCifSource::new(
        Tile::from_structure(&structure(Sample::Tiny), Seed::new(2)),
        2,
    );
    assert_ne!(drain(source(2)), drain(other));
}

#[test]
fn reading_one_byte_at_a_time_gives_the_same_bytes() {
    let whole = drain(source(2));

    let mut stream = source(2);
    let mut piecewise = Vec::new();
    let mut byte = [0_u8; 1];
    loop {
        match stream.read(&mut byte) {
            Ok(0) => break,
            Ok(_) => piecewise.push(byte[0]),
            Err(error) => panic!("byte-at-a-time read failed: {error}"),
        }
    }
    assert_eq!(String::from_utf8_lossy(&piecewise), whole);
}

#[test]
fn an_empty_target_buffer_reads_nothing() {
    let mut stream = source(1);
    match stream.read(&mut []) {
        Ok(count) => assert_eq!(count, 0),
        Err(error) => panic!("empty read failed: {error}"),
    }
}

#[test]
fn requesting_zero_copies_emits_no_coordinate_row() {
    let text = drain(source(0));
    assert_eq!(text.lines().filter(|l| l.starts_with("ATOM ")).count(), 0);
}

#[test]
fn the_stream_stays_bounded_however_long_it_runs() {
    let tile = Tile::from_structure(&structure(Sample::Tiny), Seed::new(1));
    let mut stream = SyntheticCifSource::new(tile, 512);
    let mut sink = [0_u8; 8192];
    let mut total = 0_usize;
    loop {
        match stream.read(&mut sink) {
            Ok(0) => break,
            Ok(count) => total += count,
            Err(error) => panic!("bulk read failed: {error}"),
        }
    }
    assert!(
        total > 1_000_000,
        "stream was too short to be a test: {total}"
    );
}

#[test]
fn the_generated_stream_parses_back_with_the_expected_atom_count() {
    use molframe_core::io::{InputBuffer, Limits, ReadOptions};

    let stream = source(3);
    let declared = stream.atom_count();

    let buffer = match InputBuffer::from_reader(stream, Limits::default()) {
        Ok(buffer) => buffer,
        Err(finding) => panic!("synthetic stream was refused: {finding:?}"),
    };
    let (structure, _findings) = match molframe_cif::read(&buffer, &ReadOptions::default()) {
        Ok(read) => read,
        Err(findings) => panic!("synthetic mmCIF did not parse: {findings:?}"),
    };

    assert_eq!(u64::from(structure.atom_count()), declared);
    assert!(structure.chain_count() > 0, "no chain was recovered");
    assert!(structure.residue_count() > 0, "no residue was recovered");
}

#[test]
fn a_tiled_structure_keeps_the_source_chain_count_per_copy() {
    use molframe_core::io::{InputBuffer, Limits, ReadOptions};

    let single = structure(Sample::Tiny);
    let tile = Tile::from_structure(&single, Seed::new(1));
    let copies = 3_u64;

    let buffer =
        match InputBuffer::from_reader(SyntheticCifSource::new(tile, copies), Limits::default()) {
            Ok(buffer) => buffer,
            Err(finding) => panic!("synthetic stream was refused: {finding:?}"),
        };
    let (tiled, _findings) = match molframe_cif::read(&buffer, &ReadOptions::default()) {
        Ok(read) => read,
        Err(findings) => panic!("synthetic mmCIF did not parse: {findings:?}"),
    };

    let Ok(copies) = usize::try_from(copies) else {
        panic!("copy count does not fit a usize")
    };
    assert_eq!(tiled.residue_count(), single.residue_count() * copies);
}

#[test]
fn every_copy_becomes_its_own_chain() {
    use molframe_core::io::{InputBuffer, Limits, ReadOptions};

    let single = structure(Sample::Tiny);
    let tile = Tile::from_structure(&single, Seed::new(1));

    let buffer = match InputBuffer::from_reader(SyntheticCifSource::new(tile, 4), Limits::default())
    {
        Ok(buffer) => buffer,
        Err(finding) => panic!("synthetic stream was refused: {finding:?}"),
    };
    let (tiled, _findings) = match molframe_cif::read(&buffer, &ReadOptions::default()) {
        Ok(read) => read,
        Err(findings) => panic!("synthetic mmCIF did not parse: {findings:?}"),
    };

    assert_eq!(tiled.chain_count(), single.chain_count() * 4);
}

#[test]
fn parsing_a_tiled_stream_raises_no_boundary_inference_warning() {
    use molframe_core::io::{InputBuffer, Limits, ReadOptions};

    let tile = Tile::from_structure(&structure(Sample::Tiny), Seed::new(1));
    let buffer = match InputBuffer::from_reader(SyntheticCifSource::new(tile, 3), Limits::default())
    {
        Ok(buffer) => buffer,
        Err(finding) => panic!("synthetic stream was refused: {finding:?}"),
    };
    let (_structure, findings) = match molframe_cif::read(&buffer, &ReadOptions::default()) {
        Ok(read) => read,
        Err(findings) => panic!("synthetic mmCIF did not parse: {findings:?}"),
    };

    assert!(
        findings.is_empty(),
        "the fixture should parse cleanly, but raised: {findings:?}"
    );
}
/// Streams `copies` copies and returns the bytes produced and seconds taken.
fn drain_timed(copies: u64) -> (u64, f64) {
    use std::time::Instant;

    let tile = Tile::from_structure(&structure(Sample::Large), Seed::new(1));
    let mut stream = SyntheticCifSource::new(tile, copies);
    let mut sink = vec![0_u8; 1 << 20];
    let mut total = 0_u64;

    let start = Instant::now();
    loop {
        match stream.read(&mut sink) {
            Ok(0) => break,
            Ok(count) => {
                let Ok(count) = u64::try_from(count) else {
                    panic!("read count does not fit a u64")
                };
                total += count;
            }
            Err(error) => panic!("bulk read failed: {error}"),
        }
    }
    (total, start.elapsed().as_secs_f64())
}

/// Reports the generation rate.
///
/// The rate matters — a generator slower than the reader it feeds would make
/// every measurement taken through it report the generator. It is reported here
/// rather than asserted because a throughput floor does not hold in an
/// unoptimised build, and the benchmark suite is where throughput is gated.
#[test]
#[ignore = "a throughput measurement, not an assertion; run it explicitly"]
fn report_generation_throughput() {
    let (bytes, seconds) = drain_timed(400);
    let gibibytes = crate::numeric::u64_to_f64(bytes) / crate::numeric::u64_to_f64(1 << 30);
    let rate = gibibytes / seconds;
    println!("generated {gibibytes:.2} GiB in {seconds:.2}s = {rate:.2} GiB/s");
    println!("projected wall time for 100 GiB: {:.0}s", 100.0 / rate);
}

#[test]
fn a_byte_target_produces_a_stream_of_about_that_size() {
    let tile = Tile::from_structure(&structure(Sample::Tiny), Seed::new(1));
    let target = 4_u64 << 20;
    let copies = SyntheticCifSource::copies_for_bytes(&tile, target);
    let produced = drain(SyntheticCifSource::new(tile, copies)).len();
    let Ok(produced) = u64::try_from(produced) else {
        panic!("produced length does not fit a u64")
    };
    assert!(
        produced >= target,
        "produced {produced} below target {target}"
    );
    assert!(
        produced < target * 2,
        "produced {produced} overshoots target {target}"
    );
}
