//! Deterministic provenance interchange formats.
//!
//! JSON and mmCIF are projections of one ordered field list. Keeping that list
//! in one place prevents the two representations from silently recording
//! different policy decisions.

use super::provenance::Provenance;
use std::fmt::Write;

impl Provenance {
    /// Serialises the complete record as deterministic JSON.
    ///
    /// Keys follow a fixed order and every policy field is named explicitly.
    /// The result is independent of map iteration order.
    #[must_use]
    pub fn to_json(&self) -> String {
        let fields = fields(self);
        let mut output = String::from("{");
        for (position, (key, value)) in fields.iter().enumerate() {
            if position != 0 {
                output.push(',');
            }
            json_string(&mut output, key);
            output.push(':');
            json_string(&mut output, value);
        }
        output.push('}');
        output
    }

    /// Serialises the complete record as a standalone mmCIF data block.
    ///
    /// The two-column category is deliberately lossless and extensible: older
    /// readers can preserve fields they do not interpret.
    #[must_use]
    pub fn to_mmcif(&self) -> String {
        let fields = fields(self);
        let mut output = String::from(
            "data_molframe_provenance\n#\nloop_\n_molframe_provenance.key\n_molframe_provenance.value\n",
        );
        for (key, value) in fields {
            cif_value(&mut output, &key);
            output.push(' ');
            cif_value(&mut output, &value);
            output.push('\n');
        }
        output.push_str("#\n");
        output
    }
}

fn fields(record: &Provenance) -> Vec<(String, String)> {
    let policy = &record.policy;
    let mut result = vec![
        ("molframe_version", record.molframe_version.to_owned()),
        ("input_source", record.input_source.to_string()),
        (
            "input_fingerprint",
            optional_display(record.input_fingerprint),
        ),
        ("policy_fingerprint", record.policy_fingerprint.to_string()),
        ("profile", optional_profile(record)),
        (
            "schema_version",
            optional_version(record.schema_version.as_ref()),
        ),
        (
            "component_version",
            optional_version(record.component_version.as_ref()),
        ),
        ("timestamp", optional_text(record.timestamp.as_deref())),
        ("policy.assembly", format!("{:?}", policy.assembly)),
        ("policy.model", format!("{:?}", policy.model)),
        ("policy.altloc", format!("{:?}", policy.altloc)),
        ("policy.identifiers", format!("{:?}", policy.identifiers)),
        (
            "policy.missing_atoms",
            format!("{:?}", policy.missing_atoms),
        ),
        ("policy.hydrogens", format!("{:?}", policy.hydrogens)),
        (
            "policy.atom_equivalence",
            format!("{:?}", policy.atom_equivalence),
        ),
        ("policy.symmetry", format!("{:?}", policy.symmetry)),
        ("policy.alignment", format!("{:?}", policy.alignment)),
        ("policy.precision", format!("{:?}", policy.precision)),
        ("policy.periodic", format!("{:?}", policy.periodic)),
        ("policy.vdw_radii", format!("{:?}", policy.vdw_radii)),
        ("policy.contact_def", format!("{:?}", policy.contact_def)),
        (
            "policy.float_tolerance.relative",
            policy.float_tolerance.relative.to_string(),
        ),
        (
            "policy.float_tolerance.absolute",
            policy.float_tolerance.absolute.to_string(),
        ),
    ]
    .into_iter()
    .map(|(key, value)| (key.to_owned(), value))
    .collect::<Vec<_>>();
    if let Some(algorithm) = &record.algorithm {
        result.push(("algorithm.name".to_owned(), algorithm.name().to_owned()));
        result.push((
            "algorithm.version".to_owned(),
            algorithm.version().to_owned(),
        ));
    }
    result.extend(
        record
            .parameters
            .iter()
            .map(|(name, value)| (format!("parameter.{name}"), value.serialised())),
    );
    result
}

fn optional_display(value: Option<impl std::fmt::Display>) -> String {
    match value {
        Some(value) => value.to_string(),
        None => ".".to_owned(),
    }
}

fn optional_profile(record: &Provenance) -> String {
    match record.profile {
        Some(profile) => profile.as_str().to_owned(),
        None => ".".to_owned(),
    }
}

fn optional_version(value: Option<&super::provenance::DictionaryVersion>) -> String {
    match value {
        Some(value) => value.as_str().to_owned(),
        None => ".".to_owned(),
    }
}

fn optional_text(value: Option<&str>) -> String {
    match value {
        Some(value) => value.to_owned(),
        None => ".".to_owned(),
    }
}

fn json_string(output: &mut String, value: &str) {
    output.push('"');
    for character in value.chars() {
        match character {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            control if control.is_control() => {
                let _ = write!(output, "\\u{:04x}", u32::from(control));
            }
            ordinary => output.push(ordinary),
        }
    }
    output.push('"');
}

fn cif_value(output: &mut String, value: &str) {
    if bare_cif(value) {
        output.push_str(value);
    } else if !value.contains('\'') && !value.contains('\n') {
        output.push('\'');
        output.push_str(value);
        output.push('\'');
    } else if !value.contains('"') && !value.contains('\n') {
        output.push('"');
        output.push_str(value);
        output.push('"');
    } else {
        output.push_str("\n;");
        output.push_str(value);
        output.push_str("\n;");
    }
}

fn bare_cif(value: &str) -> bool {
    !value.is_empty()
        && value != "."
        && value != "?"
        && !value.chars().any(char::is_whitespace)
        && !value.starts_with(['#', '$', ';', '_', '\'', '"'])
        && !value.eq_ignore_ascii_case("loop_")
        && !value.to_ascii_lowercase().starts_with("data_")
        && !value.to_ascii_lowercase().starts_with("save_")
}

#[cfg(test)]
#[path = "serialise_tests.rs"]
mod tests;
