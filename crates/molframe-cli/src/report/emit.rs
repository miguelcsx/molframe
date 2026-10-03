//! Row emission shared by every command that prints a table.

use super::{Context, RowWriter};
use crate::exit::Exit;

/// Writes one row per item, rendering each with `render`.
pub(crate) fn emit_rows<T>(
    context: Context,
    header: &[&str],
    rows: impl IntoIterator<Item = T>,
    render: impl Fn(T) -> Vec<String>,
) -> Exit {
    let mut output = match RowWriter::new(context, header) {
        Ok(output) => output,
        Err(error) => return output_error(&error),
    };
    for item in rows {
        let row = render(item);
        if let Err(error) = output.row(row.iter().map(String::as_str)) {
            return output_error(&error);
        }
    }
    finish_rows(output, None)
}

pub(crate) fn finish_rows(output: RowWriter, error: Option<&std::io::Error>) -> Exit {
    if let Some(error) = error {
        return output_error(error);
    }
    match output.finish() {
        Ok(()) => Exit::Success,
        Err(error) => output_error(&error),
    }
}

pub(crate) fn output_error(error: &std::io::Error) -> Exit {
    eprintln!("could not write result: {error}");
    Exit::Consistency
}
