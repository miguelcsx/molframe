//! Grace/xmgrace time-series tables, as simulation tools write them.
//!
//! An energy log, an RMSD trace and a per-residue distance plot all arrive in
//! this format: comment lines beginning `#`, plot-metadata lines beginning `@`,
//! and then one row of numbers per record with the abscissa first. Nothing here
//! interprets the columns — the file says what they are through its legends,
//! and the caller decides what to do with them.
//!
//! Reading is one pass over the text, `O(values)`. Columns are stored one
//! contiguous vector each rather than a vector of rows: a consumer wants a
//! whole series at a time, and a column layout also lets a single series be
//! handed onward as a slice without a gather.

use std::collections::BTreeMap;

/// A parsed time-series table.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Series {
    /// Plot title, when the file declared one.
    pub title: Option<Box<str>>,
    /// Abscissa axis label, typically a time unit.
    pub abscissa_label: Option<Box<str>>,
    /// Ordinate axis label.
    pub ordinate_label: Option<Box<str>>,
    /// Legend for each ordinate column, in column order; absent entries are
    /// empty.
    pub legends: Vec<Box<str>>,
    /// Abscissa value of every record, usually simulation time.
    pub abscissa: Vec<f64>,
    /// One contiguous vector per ordinate column.
    pub columns: Vec<Vec<f64>>,
}

impl Series {
    /// Records in the table.
    #[must_use]
    pub fn len(&self) -> usize {
        self.abscissa.len()
    }

    /// Whether the table holds no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.abscissa.is_empty()
    }

    /// The column whose legend matches exactly, if any.
    ///
    /// Legends are how a simulation tool names what it wrote, so matching on
    /// one is how a caller asks for "Potential" or "RMSD" without hard-coding
    /// a column index that changes with the tool's version. `O(columns)`.
    #[must_use]
    pub fn column_by_legend(&self, legend: &str) -> Option<&[f64]> {
        let index = self
            .legends
            .iter()
            .position(|candidate| candidate.as_ref() == legend)?;
        self.columns.get(index).map(Vec::as_slice)
    }

    /// Linearly interpolates one column at an arbitrary abscissa.
    ///
    /// A trajectory and its energy log rarely share a sample rate, and a
    /// viewer showing both at one timestamp needs a value between records
    /// rather than the nearest one. Values outside the recorded range clamp to
    /// the first or last record. `O(log records)`, since the abscissa of these
    /// files is monotonic by construction.
    #[must_use]
    pub fn sample(&self, column: usize, abscissa: f64) -> Option<f64> {
        let values = self.columns.get(column)?;
        if values.is_empty() || values.len() != self.abscissa.len() {
            return None;
        }
        let position = self
            .abscissa
            .partition_point(|candidate| *candidate < abscissa);
        if position == 0 {
            return values.first().copied();
        }
        if position >= values.len() {
            return values.last().copied();
        }
        let (before, after) = (self.abscissa[position - 1], self.abscissa[position]);
        let span = after - before;
        if span <= 0.0 {
            return values.get(position).copied();
        }
        let weight = (abscissa - before) / span;
        Some(values[position - 1] + weight * (values[position] - values[position - 1]))
    }
}

/// Why a time-series table could not be read.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum SeriesError {
    /// The text is not UTF-8.
    #[error("time-series table is not valid UTF-8")]
    NotUtf8,
    /// A data row held a token that is not a number.
    #[error("time-series row {row} column {column} is not a number")]
    InvalidNumber {
        /// Zero-based data row.
        row: usize,
        /// Zero-based column within the row.
        column: usize,
    },
    /// A data row held a different number of columns than the first row.
    #[error("time-series row {row} has {found} columns, expected {expected}")]
    RaggedRow {
        /// Zero-based data row.
        row: usize,
        /// Columns found on this row.
        found: usize,
        /// Columns the first row established.
        expected: usize,
    },
    /// The host could not reserve the table.
    #[error("time-series table allocation exceeds host resources")]
    ResourceLimit,
}

/// Reads a Grace time-series table.
///
/// # Errors
///
/// Returns [`SeriesError`] for non-UTF-8 text, a row whose width disagrees with
/// the first, a token that is not a number, or a table the host cannot hold.
pub fn read_xvg(bytes: &[u8]) -> Result<Series, SeriesError> {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return Err(SeriesError::NotUtf8);
    };
    let mut series = Series::default();
    let mut legends: BTreeMap<usize, Box<str>> = BTreeMap::new();
    let mut width = None;
    let mut row = 0;

    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('@') {
            read_metadata(rest.trim(), &mut series, &mut legends);
            continue;
        }
        // A dataset separator; the format allows several sets per file and this
        // reader keeps them as one continuous table.
        if trimmed.starts_with('&') {
            continue;
        }

        let mut values = trimmed.split_ascii_whitespace();
        let Some(first) = values.next() else {
            continue;
        };
        let Ok(abscissa) = first.parse::<f64>() else {
            return Err(SeriesError::InvalidNumber { row, column: 0 });
        };

        let mut ordinates = 0;
        for (index, token) in values.enumerate() {
            let Ok(value) = token.parse::<f64>() else {
                return Err(SeriesError::InvalidNumber {
                    row,
                    column: index + 1,
                });
            };
            if series.columns.len() <= index {
                // The first data row fixes the table's width; a later row with
                // an extra column is a ragged file, not a wider table.
                if let Some(expected) = width {
                    return Err(SeriesError::RaggedRow {
                        row,
                        found: index + 2,
                        expected,
                    });
                }
                series.columns.push(Vec::new());
            }
            series.columns[index].push(value);
            ordinates = index + 1;
        }

        match width {
            None => width = Some(ordinates + 1),
            Some(expected) if expected != ordinates + 1 => {
                return Err(SeriesError::RaggedRow {
                    row,
                    found: ordinates + 1,
                    expected,
                });
            }
            Some(_) => {}
        }
        series.abscissa.push(abscissa);
        row += 1;
    }

    series.legends = (0..series.columns.len())
        .map(|index| match legends.remove(&index) {
            Some(legend) => legend,
            None => Box::from(""),
        })
        .collect();
    Ok(series)
}

/// Interprets the `@` metadata lines this reader cares about, ignoring the
/// rest of Grace's plot-styling vocabulary.
fn read_metadata(rest: &str, series: &mut Series, legends: &mut BTreeMap<usize, Box<str>>) {
    if let Some(value) = quoted_after(rest, "title") {
        series.title = Some(value);
        return;
    }
    if let Some(value) = quoted_after(rest, "xaxis  label") {
        series.abscissa_label = Some(value);
        return;
    }
    if let Some(value) = quoted_after(rest, "xaxis label") {
        series.abscissa_label = Some(value);
        return;
    }
    if let Some(value) = quoted_after(rest, "yaxis  label") {
        series.ordinate_label = Some(value);
        return;
    }
    if let Some(value) = quoted_after(rest, "yaxis label") {
        series.ordinate_label = Some(value);
        return;
    }
    // Series legends arrive as `s0 legend "Potential"`.
    let Some(after_s) = rest.strip_prefix('s') else {
        return;
    };
    let digits: String = after_s.chars().take_while(char::is_ascii_digit).collect();
    let Ok(index) = digits.parse::<usize>() else {
        return;
    };
    if let Some(value) = quoted_after(&after_s[digits.len()..], "legend") {
        legends.insert(index, value);
    }
}

/// The quoted string following a keyword, when the line begins with it.
fn quoted_after(rest: &str, keyword: &str) -> Option<Box<str>> {
    let trimmed = rest.trim_start();
    let after = trimmed.strip_prefix(keyword)?;
    let start = after.find('"')?;
    let remainder = &after[start + 1..];
    let end = remainder.find('"')?;
    Some(Box::from(&remainder[..end]))
}

#[cfg(test)]
#[path = "xvg_tests.rs"]
mod tests;
