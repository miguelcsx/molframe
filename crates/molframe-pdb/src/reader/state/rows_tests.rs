//! What a coordinate-looking row that will not parse is reported as.

use super::starts_with_coordinate_keyword;

/// A well-formed row names a record in columns 1 to 6 and is not a layout fault.
#[test]
fn a_well_formed_coordinate_row_is_not_malformed() {
    assert!(!starts_with_coordinate_keyword(
        "ATOM      1  N   ALA A   1      27.340  24.430   2.614  1.00 10.00           N"
    ));
    assert!(!starts_with_coordinate_keyword(
        "HETATM    1  O   HOH A 100      27.340  24.430   2.614  1.00 10.00           O"
    ));
}

/// A serial that outgrew its five columns shifts a digit into the record name.
///
/// `ATOM` is four characters, so a six-digit serial reaches column 6 and the name
/// reads `ATOM 1`. `HETATM` already fills six columns, so it cannot shift into its
/// own name and instead pushes the rest of the row right.
#[test]
fn a_serial_that_outgrew_its_field_is_recognised_as_a_coordinate_row() {
    assert!(starts_with_coordinate_keyword(
        "ATOM 100000  N   ALA A   1      27.340  24.430   2.614  1.00 10.00           N"
    ));
    // An `ATOM` name with two digits shifted in.
    assert!(starts_with_coordinate_keyword(
        "ATOM 1000000  N   ALA A   1      27.340  24.430   2.614  1.00 10.00          N"
    ));
}

/// A line that is not a coordinate row at all stays out of the count.
#[test]
fn other_records_and_prose_are_not_coordinate_rows() {
    for line in [
        "HEADER    HYDROLASE                               01-JAN-00   1ABC",
        "TER       7      ASN A   2",
        "REMARK   2 RESOLUTION.    1.50 ANGSTROMS.",
        "",
        "ATOM",
    ] {
        assert!(
            !starts_with_coordinate_keyword(line),
            "`{line}` is not a coordinate row"
        );
    }
}

/// A coordinate record whose serial outgrew its five columns is reported as a
/// field that would not parse, not as a file with no coordinates.
///
/// The record name is columns 1 to 6, so `ATOM 100000` reads `ATOM 1` and is not
/// an `ATOM` row at all; the format requires the hybrid-36 spelling `A0000` for
/// that serial. Saying "the file contains no coordinate records" about a file
/// that plainly holds one leaves a caller nothing to act on, so the reader counts
/// the rows whose keyword is present but whose fields do not line up.
#[test]
fn a_serial_that_outgrew_its_field_is_reported_as_an_unparseable_row() {
    use molframe_core::diagnostic::Code;
    use molframe_core::io::{InputBuffer, ReadOptions};

    let input = InputBuffer::from_bytes(
        "ATOM 100000  N   ALA A   1      27.340  24.430   2.614  1.00 10.00           N\n\
END\n"
            .as_bytes()
            .to_vec(),
    );
    let Err(findings) = crate::read(&input, &ReadOptions::new()) else {
        panic!("a serial outside the field layout cannot be read as an atom row");
    };
    assert!(
        findings.iter().any(|finding| {
            finding.code() == Code::E1201 && finding.message().contains("fields do not parse")
        }),
        "the report must name the field layout, not claim the file is empty: {findings:?}"
    );
    // And the hybrid-36 spelling of the same serial is a valid row.
    let input = InputBuffer::from_bytes(
        "ATOM  A0000  N   ALA A   1      27.340  24.430   2.614  1.00 10.00           N\n\
END\n"
            .as_bytes()
            .to_vec(),
    );
    let (structure, _) = crate::read(&input, &ReadOptions::new())
        .unwrap_or_else(|findings| panic!("hybrid-36 serial: {findings:?}"));
    assert_eq!(structure.atom_count(), 1);
}

/// A file that genuinely holds no coordinates still says so.
#[test]
fn a_file_without_coordinate_records_reports_exactly_that() {
    use molframe_core::diagnostic::Code;
    use molframe_core::io::{InputBuffer, ReadOptions};

    let input = InputBuffer::from_bytes(b"HEADER    NOTHING\nEND\n".to_vec());
    let Err(findings) = crate::read(&input, &ReadOptions::new()) else {
        panic!("no coordinates is a failure");
    };
    assert!(
        findings.iter().any(|finding| finding.code() == Code::E1001
            && finding.message().contains("no coordinate records")),
        "{findings:?}"
    );
    assert!(
        !findings.iter().any(|finding| finding.code() == Code::E1201),
        "a file with no coordinate rows at all must not be reported as a layout fault: {findings:?}"
    );
}
