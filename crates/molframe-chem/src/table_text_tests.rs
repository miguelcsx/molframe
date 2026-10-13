use super::*;

struct Collect(Vec<(usize, String)>, bool);

impl TableLines for Collect {
    fn header(&mut self, _line: usize, key: &str, _value: &str) -> Result<bool, TableTextError> {
        self.1 = key == "base";
        Ok(self.1)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        self.0.push((line, text.to_owned()));
        Ok(())
    }
}

const HEAD: &str = "# SPDX-License-Identifier: CC0-1.0\nformat: 1\nname: n\nversion: 1\nsource: s\nlicence: CC0-1.0\n";

fn read(text: &str) -> (Result<TableIdentity, TableTextError>, Collect) {
    let mut sink = Collect(Vec::new(), false);
    (read_table("test table", text, &mut sink), sink)
}

#[test]
fn headers_become_the_identity_and_the_rest_are_rows() {
    let (identity, sink) = read(&format!("{HEAD}base: x\nA 1\n"));
    assert_eq!(
        identity,
        Ok(TableIdentity::supplied("n", "1", "s", "CC0-1.0"))
    );
    assert_eq!(sink.0, vec![(8, "A 1".to_owned())]);
    assert!(sink.1);
}

#[test]
fn the_format_line_is_optional_but_only_version_one_is_understood() {
    assert!(
        read("name: n\nversion: 1\nsource: s\nlicence: l\n")
            .0
            .is_ok()
    );
    let (bad, _) = read(&HEAD.replace("format: 1", "format: 2"));
    assert_eq!(bad.map_err(|error| error.line()), Err(2));
}

#[test]
fn a_missing_repeated_or_unknown_header_is_refused() {
    let (missing, _) = read("name: n\nversion: 1\nsource: s\n");
    assert_eq!(missing.map_err(|error| error.line()), Err(0));
    let (twice, _) = read(&format!("{HEAD}name: m\n"));
    assert_eq!(twice.map_err(|error| error.line()), Err(7));
    let (unknown, _) = read(&format!("{HEAD}colour: red\n"));
    assert_eq!(unknown.map_err(|error| error.line()), Err(7));
}

#[test]
fn a_header_without_a_value_is_refused() {
    let (empty, _) = read("name:\n");
    assert_eq!(empty.map_err(|error| error.line()), Err(1));
}

#[test]
fn bytes_that_are_not_text_are_refused() {
    assert!(text_of("t", &[0xff, 0xfe]).is_err());
    assert_eq!(text_of("t", b"ok"), Ok("ok"));
}
