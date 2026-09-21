//! Incremental row output with atomic file publication.

use super::{Context, OutputKind, delimited_field, escape};
use std::fmt::Write as _;
use std::io::{self, Write as _};
use std::path::{Path, PathBuf};

/// A bounded row sink that never retains prior rows.
pub struct RowWriter {
    context: Context,
    header: Vec<String>,
    destination: Option<Destination>,
    first: bool,
    buffer: String,
}

enum Destination {
    Stdout(std::io::Stdout),
    File {
        temporary: tempfile::NamedTempFile,
        path: PathBuf,
    },
}

impl RowWriter {
    /// Opens a row stream and emits its format-specific prefix or header.
    pub fn new(context: Context, header: &[&str]) -> io::Result<Self> {
        let destination = match context.output {
            Some(path) => Destination::File {
                temporary: temporary_for(path)?,
                path: path.to_owned(),
            },
            None => Destination::Stdout(std::io::stdout()),
        };
        let mut writer = Self {
            context,
            header: header.iter().map(|value| (*value).to_owned()).collect(),
            destination: Some(destination),
            first: true,
            buffer: String::new(),
        };
        writer.begin()?;
        Ok(writer)
    }

    /// Emits one row without retaining it.
    pub fn row<'a>(&mut self, values: impl IntoIterator<Item = &'a str>) -> io::Result<()> {
        let mut row = std::mem::take(&mut self.buffer);
        row.clear();
        let mut width = 0;
        match self.context.format {
            OutputKind::Json | OutputKind::JsonLines => {
                if self.context.format == OutputKind::Json && !self.first {
                    row.push(',');
                }
                row.push('{');
                for (position, value) in values.into_iter().enumerate() {
                    let Some(key) = self.header.get(position) else {
                        self.buffer = row;
                        return Err(width_error());
                    };
                    if position != 0 {
                        row.push(',');
                    }
                    let _ = write!(row, "\"{}\":\"{}\"", escape(key), escape(value));
                    width += 1;
                }
                row.push('}');
                if self.context.format == OutputKind::JsonLines {
                    row.push('\n');
                }
            }
            OutputKind::Csv | OutputKind::Tsv | OutputKind::Text => {
                let delimiter = self.context.table_delimiter();
                for (position, value) in values.into_iter().enumerate() {
                    if position != 0 {
                        row.push(delimiter);
                    }
                    delimited_field(&mut row, value, delimiter);
                    width += 1;
                }
                row.push('\n');
            }
            OutputKind::Arrow | OutputKind::Parquet => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "row sink cannot encode Arrow or Parquet",
                ));
            }
        }
        if width != self.header.len() {
            self.buffer = row;
            return Err(width_error());
        }
        self.first = false;
        let result = self.write_all(row.as_bytes());
        row.clear();
        self.buffer = row;
        result
    }

    /// Emits one pre-encoded JSON object, preserving numeric field types.
    pub fn json_record(&mut self, record: &str) -> io::Result<()> {
        if !self.context.is_json() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "JSON record needs JSON output",
            ));
        }
        if self.context.format == OutputKind::Json && !self.first {
            self.write_all(b",")?;
        }
        self.write_all(record.as_bytes())?;
        if self.context.format == OutputKind::JsonLines {
            self.write_all(b"\n")?;
        }
        self.first = false;
        Ok(())
    }

    /// Completes the stream and atomically publishes file output.
    pub fn finish(mut self) -> io::Result<()> {
        let mut suffix = String::new();
        match self.context.format {
            OutputKind::Json => {
                let _ = writeln!(
                    suffix,
                    "],\"_provenance\":{}}}",
                    self.context.provenance_json()
                );
            }
            OutputKind::JsonLines => {
                let _ = writeln!(
                    suffix,
                    "{{\"_provenance\":{}}}",
                    self.context.provenance_json()
                );
            }
            OutputKind::Csv | OutputKind::Tsv => {
                let _ = writeln!(
                    suffix,
                    "# molframe-provenance={}",
                    self.context.provenance_json()
                );
            }
            OutputKind::Text => {}
            OutputKind::Arrow | OutputKind::Parquet => unreachable!(),
        }
        self.write_all(suffix.as_bytes())?;
        self.flush()?;
        self.publish()?;
        if let Some(path) = self.context.provenance {
            atomic_write(
                path,
                format!("{}\n", self.context.provenance_json()).as_bytes(),
            )?;
        }
        Ok(())
    }

    fn begin(&mut self) -> io::Result<()> {
        match self.context.format {
            OutputKind::Json => self.write_all(b"{\"result\":["),
            OutputKind::JsonLines => Ok(()),
            OutputKind::Csv | OutputKind::Tsv | OutputKind::Text => {
                let delimiter = self.context.table_delimiter();
                let mut line = String::new();
                append_delimited(&mut line, self.header.iter().map(String::as_str), delimiter);
                self.write_all(line.as_bytes())
            }
            OutputKind::Arrow | OutputKind::Parquet => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "row sink cannot encode Arrow or Parquet",
            )),
        }
    }

    fn write_all(&mut self, bytes: &[u8]) -> io::Result<()> {
        match self.destination.as_mut() {
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "row sink is closed",
            )),
            Some(Destination::Stdout(output)) => output.write_all(bytes),
            Some(Destination::File { temporary, .. }) => temporary.write_all(bytes),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self.destination.as_mut() {
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "row sink is closed",
            )),
            Some(Destination::Stdout(output)) => output.flush(),
            Some(Destination::File { temporary, .. }) => temporary.flush(),
        }
    }

    fn publish(&mut self) -> io::Result<()> {
        match self.destination.take() {
            Some(Destination::Stdout(output)) => {
                self.destination = Some(Destination::Stdout(output));
                Ok(())
            }
            Some(Destination::File { temporary, path }) => {
                temporary.persist(path).map_err(|error| error.error)?;
                Ok(())
            }
            None => Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "row sink is closed",
            )),
        }
    }
}

fn append_delimited<'a>(
    output: &mut String,
    values: impl IntoIterator<Item = &'a str>,
    delimiter: char,
) {
    for (position, value) in values.into_iter().enumerate() {
        if position != 0 {
            output.push(delimiter);
        }
        delimited_field(output, value, delimiter);
    }
    output.push('\n');
}

fn width_error() -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        "row width does not match output header",
    )
}

fn temporary_for(path: &Path) -> io::Result<tempfile::NamedTempFile> {
    match path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        Some(parent) => tempfile::NamedTempFile::new_in(parent),
        None => tempfile::NamedTempFile::new_in("."),
    }
}

fn atomic_write(path: &Path, bytes: &[u8]) -> io::Result<()> {
    let mut temporary = temporary_for(path)?;
    temporary.write_all(bytes)?;
    temporary.flush()?;
    temporary.persist(path).map_err(|error| error.error)?;
    Ok(())
}
