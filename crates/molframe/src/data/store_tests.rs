use super::*;
use std::io::{Read as _, Write as _};
use std::net::TcpListener;
use std::time::Duration;

const BODY: &str = "Bond\tResidue\tMean\tStdDev\nCA-CB\tALA\t1.52\t0.021\n-\n";

fn options() -> DownloadOptions {
    DownloadOptions {
        max_bytes: 1 << 16,
        timeout: Duration::from_secs(5),
        redirect_limit: 0,
    }
}

fn sha(bytes: &[u8]) -> String {
    use sha2::{Digest as _, Sha256};
    use std::fmt::Write as _;
    Sha256::digest(bytes)
        .iter()
        .fold(String::new(), |mut out, byte| {
            let _ = write!(out, "{byte:02x}");
            out
        })
}

fn serve(body: &'static str) -> String {
    let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
    let address = listener.local_addr().expect("address");
    std::thread::spawn(move || {
        if let Ok((mut stream, _)) = listener.accept() {
            let mut request = [0_u8; 1024];
            let _ = stream.read(&mut request);
            let head = format!(
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = stream.write_all(head.as_bytes());
            let _ = stream.write_all(body.as_bytes());
        }
    });
    format!("http://{address}/props.txt")
}

fn entry(url: String, append: Vec<PathBuf>) -> Entry {
    Entry {
        name: "geo".to_owned(),
        version: "1".to_owned(),
        url,
        licence: "MIT".to_owned(),
        sha256: sha(BODY.as_bytes()),
        extension: "geometry".to_owned(),
        source: Some("synthetic".to_owned()),
        converter: Some("stereo_chemical_props".to_owned()),
        append,
    }
}

#[test]
fn fetch_converts_appends_stores_and_verifies() {
    let directory = tempfile::tempdir().expect("directory");
    let extra = directory.path().join("extra.rows");
    std::fs::write(&extra, "# extra\nlink\t*\tC-N\t1.329\t0.014\n").expect("write");
    let manifest = Manifest::from_entries(
        directory.path().join("manifest.toml"),
        vec![entry(serve(BODY), vec![PathBuf::from("extra.rows")])],
    );
    let cache = directory.path().join("cache");
    assert_eq!(list(&manifest, &cache)[0].status, Status::Missing);
    let stored = fetch(&manifest, "geo", &cache, false, options()).expect("fetched");
    assert_eq!(stored.rows, 2);
    let text = std::fs::read_to_string(&stored.path).expect("read");
    assert!(text.contains("bond\tALA\tCA-CB\t1.52\t0.021\n") && text.contains("link\t*\tC-N"));
    assert_eq!(
        verify(&manifest, &cache, None).expect("verify")[0].status,
        Status::Current
    );
    assert!(matches!(
        fetch(&manifest, "geo", &cache, false, options()),
        Err(DataError::Exists { .. })
    ));
}

#[test]
fn a_changed_pin_makes_the_stored_table_stale_and_garbage_invalid() {
    let directory = tempfile::tempdir().expect("directory");
    let manifest = Manifest::from_entries(
        directory.path().join("m.toml"),
        vec![entry("http://127.0.0.1:1/x".to_owned(), Vec::new())],
    );
    let path = directory.path().join("geo.geometry");
    let rows = ["bond\tALA\tA-B\t1\t0.1".to_owned()];
    std::fs::write(&path, table_text(&manifest.entries()[0], &rows)).expect("write");
    assert_eq!(
        verify(&manifest, directory.path(), Some("geo")).expect("verify")[0].status,
        Status::Current
    );
    let mut moved = manifest.entries()[0].clone();
    moved.sha256 = "1".repeat(64);
    let repinned = Manifest::from_entries(directory.path().join("m.toml"), vec![moved]);
    assert!(matches!(
        verify(&repinned, directory.path(), None).expect("verify")[0].status,
        Status::Stale(_)
    ));
    std::fs::write(&path, "name: geo\n").expect("write");
    assert!(matches!(
        verify(&manifest, directory.path(), None).expect("verify")[0].status,
        Status::Invalid(_)
    ));
}

#[test]
fn unknown_entries_converters_and_append_files_fail_before_any_download() {
    let directory = tempfile::tempdir().expect("directory");
    let mut bad = entry(
        "http://127.0.0.1:1/x".to_owned(),
        vec![PathBuf::from("none.rows")],
    );
    let manifest = Manifest::from_entries(directory.path().join("m.toml"), vec![bad.clone()]);
    assert!(matches!(
        fetch(&manifest, "other", directory.path(), false, options()),
        Err(DataError::NoEntry { .. })
    ));
    assert!(matches!(
        fetch(&manifest, "geo", directory.path(), false, options()),
        Err(DataError::Read { .. })
    ));
    bad.converter = Some("nope".to_owned());
    let manifest = Manifest::from_entries(directory.path().join("m.toml"), vec![bad]);
    assert!(matches!(
        fetch(&manifest, "geo", directory.path(), false, options()),
        Err(DataError::Convert { .. })
    ));
}
