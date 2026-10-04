use super::certificate;
use crate::{PolicyDimension, PolicySpace, ScalarError, ScalarMode, audit_analyses};
use molframe_core::contract::{
    AlgorithmId, Analysis, AnalysisPolicy, ContentDigest, Coverage, HydrogenPolicy, PolicyField,
    SourceRef,
};
use std::convert::Infallible;

fn plan() -> crate::AuditPlan {
    let dimension =
        PolicyDimension::hydrogens([HydrogenPolicy::ExplicitOnly, HydrogenPolicy::Exclude])
            .justified(
                "prepared files add hydrogens that deposited ones lack",
                "PDBbind v2020 protein files",
            );
    match PolicySpace::new(AnalysisPolicy::default())
        .vary(dimension)
        .plan()
    {
        Ok(plan) => plan,
        Err(error) => panic!("plan failed: {error}"),
    }
}

fn run(digest: Option<ContentDigest>) -> String {
    let metric = ScalarError::new(|value: &f64| *value, ScalarMode::Absolute);
    let analyse = |policy: &AnalysisPolicy| -> Result<Analysis<f64>, Infallible> {
        let value = if matches!(policy.hydrogens, HydrogenPolicy::Exclude) {
            1.0
        } else {
            3.0
        };
        let mut result = Analysis::complete(value, Coverage::complete(1), policy);
        result.provenance = result
            .provenance
            .with_policy_reads(&[PolicyField::Hydrogens])
            .with_source(SourceRef::path("1abc.cif"))
            .with_algorithm(AlgorithmId::new("atom-contacts", "1"))
            .with_estimand("the atom pairs within the cutoff");
        if let Some(digest) = digest {
            result.provenance = result.provenance.with_input_digest(digest);
        }
        Ok(result)
    };
    let plan = plan();
    let Ok(audit) = audit_analyses(&plan, analyse, &metric) else {
        panic!("a plan over a field the analysis reads");
    };
    certificate(&plan, &audit)
}

/// Whether braces and brackets balance outside strings, which every well-formed
/// document does; the Python suite parses the document properly.
fn balanced(document: &str) -> bool {
    let (mut depth, mut inside, mut escaped) = (0_i64, false, false);
    for character in document.chars() {
        if inside {
            match (escaped, character) {
                (true, _) => escaped = false,
                (false, '\\') => escaped = true,
                (false, '"') => inside = false,
                _ => {}
            }
            continue;
        }
        match character {
            '"' => inside = true,
            '{' | '[' => depth += 1,
            '}' | ']' => depth -= 1,
            _ => {}
        }
        if depth < 0 {
            return false;
        }
    }
    depth == 0 && !inside
}

#[test]
fn the_certificate_names_every_universe_decision_and_finding() {
    let digest = ContentDigest::of(b"data_1abc");
    let document = run(Some(digest));
    assert!(balanced(&document));
    assert!(document.starts_with("{\"@context\":\"https://w3id.org/ro/crate/1.1/context\""));
    assert!(document.contains("\"@id\":\"#universe-0\""));
    assert!(document.contains("\"@id\":\"#universe-1\""));
    assert!(document.contains("\"@id\":\"#decision-hydrogens\""));
    assert!(document.contains("prepared files add hydrogens"));
    assert!(document.contains("PDBbind v2020 protein files"));
    assert!(document.contains("\"uncertainty_class\""));
    assert!(document.contains(&digest.to_string()));
    assert!(document.contains("the atom pairs within the cutoff"));
    assert!(document.contains("\"name\":\"main_effect.hydrogens.share\""));
    assert!(document.contains("does not say which choice is correct"));
}

#[test]
fn the_same_audit_of_the_same_bytes_writes_the_same_document() {
    let digest = ContentDigest::of(b"data_1abc");
    assert_eq!(run(Some(digest)), run(Some(digest)));
}

#[test]
fn an_input_without_a_digest_is_counted_not_vouched_for() {
    let document = run(None);
    assert!(document.contains("\"name\":\"inputs_without_sha256\",\"value\":1"));
    assert!(!document.contains("\"name\":\"sha256\""));
}
