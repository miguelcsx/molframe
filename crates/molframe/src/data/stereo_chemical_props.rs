//! Converter for a stereochemical parameter file.
//!
//! The file has a `Bond` section and an `Angle` section, each introduced by a
//! header line and closed by a line holding `-`, with rows of
//! `atoms residue mean sd`. A later section of non-bonded distances is skipped.
//! Rows become `bond|angle <TAB> residue <TAB> atoms <TAB> mean <TAB> sd`, in
//! file order. Anything else inside a section is refused with its line number.
//! Cost: one pass, one allocation per row.

#[derive(Clone, Copy, PartialEq, Eq)]
enum Section {
    Outside,
    Bond,
    Angle,
    Skipped,
}

/// The geometry rows of `text`.
pub(super) fn rows(text: &str) -> Result<Vec<String>, String> {
    let mut section = Section::Outside;
    let mut out = Vec::new();
    for (index, raw) in text.lines().enumerate() {
        let number = index + 1;
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line == "-" {
            section = Section::Outside;
            continue;
        }
        let Some(first) = line.split_whitespace().next() else {
            continue;
        };
        if section == Section::Outside || first == "Bond" || first == "Angle" {
            section = match first {
                "Bond" => Section::Bond,
                "Angle" => Section::Angle,
                _ => Section::Skipped,
            };
            continue;
        }
        let (kind, atoms) = match section {
            Section::Bond => ("bond", 2),
            Section::Angle => ("angle", 3),
            Section::Outside | Section::Skipped => continue,
        };
        out.push(row(number, line, kind, atoms)?);
    }
    if out.is_empty() {
        return Err("no bond or angle rows found".to_owned());
    }
    Ok(out)
}

fn row(number: usize, line: &str, kind: &str, atoms: usize) -> Result<String, String> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    let [names, residue, mean, sd] = fields.as_slice() else {
        return Err(format!(
            "line {number}: expected `atoms residue mean sd`, found `{line}`"
        ));
    };
    if names.split('-').count() != atoms || names.split('-').any(str::is_empty) {
        return Err(format!(
            "line {number}: `{names}` is not {atoms} atom names joined by `-`"
        ));
    }
    for (label, value) in [("mean", mean), ("sd", sd)] {
        if !value.parse::<f64>().is_ok_and(f64::is_finite) {
            return Err(format!("line {number}: {label} `{value}` is not a number"));
        }
    }
    Ok(format!("{kind}\t{residue}\t{names}\t{mean}\t{sd}"))
}

#[cfg(test)]
#[path = "stereo_chemical_props_tests.rs"]
mod tests;
