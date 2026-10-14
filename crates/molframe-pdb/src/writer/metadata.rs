//! Header records the writer generates from entry metadata.

use super::PdbOptions;
use molframe_core::structure::Structure;
use std::fmt;

const FIRST_TITLE_WIDTH: usize = 70;
const NEXT_TITLE_WIDTH: usize = 69;

/// Write HEADER, TITLE, EXPDTA, resolution, cell, and caller-supplied records.
pub(super) fn write_generated(
    out: &mut impl fmt::Write,
    structure: &Structure,
    options: &PdbOptions,
) {
    let data = structure.data();
    if let Some(id) = &data.entry.id {
        let id: String = id.chars().take(4).collect();
        let _ = writeln!(out, "HEADER    {:<40}{:<9}   {:>4}", "", "", id);
    }
    if let Some(title) = &data.entry.title {
        write_title(out, title);
    }
    if let Some(method) = &data.entry.method {
        let _ = writeln!(out, "EXPDTA    {method}");
    }
    if let Some(resolution) = data.entry.resolution {
        let _ = writeln!(out, "REMARK   2");
        let _ = writeln!(out, "REMARK   2 RESOLUTION. {resolution:>7.2} ANGSTROMS.");
    }
    if let Some(cell) = data.cell {
        let symbol = match data
            .entry
            .space_group
            .as_deref()
            .or_else(|| options.space_group_fallback())
        {
            Some(symbol) => symbol,
            None => "P 1",
        };
        let _ = writeln!(
            out,
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:<11}",
            cell.lengths[0],
            cell.lengths[1],
            cell.lengths[2],
            cell.angles[0],
            cell.angles[1],
            cell.angles[2],
            symbol,
        );
        super::connections::write_scale(out, &cell);
    }
    if let Some(extra) = options.extra_records() {
        for line in extra.lines() {
            let _ = writeln!(out, "{line}");
        }
    }
}

/// Wrap a title into continuation lines at word boundaries.
fn write_title(out: &mut impl fmt::Write, title: &str) {
    let mut line = String::new();
    let mut number = 1usize;
    let flush = |out: &mut dyn fmt::Write, line: &mut String, number: &mut usize| {
        if line.is_empty() {
            return;
        }
        if *number == 1 {
            let _ = writeln!(out, "TITLE     {line}");
        } else {
            let _ = writeln!(out, "TITLE   {:>2} {line}", *number);
        }
        *number += 1;
        line.clear();
    };
    for word in title.split_whitespace() {
        let mut rest = word;
        loop {
            let width = if number == 1 {
                FIRST_TITLE_WIDTH
            } else {
                NEXT_TITLE_WIDTH
            };
            let used = line.chars().count();
            let separator = usize::from(used > 0);
            let remaining = rest.chars().count();
            if used + separator + remaining <= width {
                if separator == 1 {
                    line.push(' ');
                }
                line.push_str(rest);
                break;
            }
            if used > 0 {
                flush(out, &mut line, &mut number);
                continue;
            }
            let split = rest
                .char_indices()
                .nth(width)
                .map_or(rest.len(), |(i, _)| i);
            line.push_str(&rest[..split]);
            rest = &rest[split..];
            flush(out, &mut line, &mut number);
        }
    }
    flush(out, &mut line, &mut number);
}

#[cfg(test)]
#[path = "metadata_tests.rs"]
mod tests;
