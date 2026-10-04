//! A semantic robustness certificate: what was asked, of which bytes, under which
//! decisions, and how far the answer moved.
//!
//! The certificate is an RO-Crate 1.1 metadata document (JSON-LD). It uses the
//! schema.org vocabulary a process-run crate uses: the software and the algorithm are
//! `SoftwareApplication`s, each universe is a `CreateAction` with the exact policy as
//! properties, the audit as a whole is an `OrganizeAction` over them, and what it
//! measured is a `Dataset` of named values. Mapped to PROV, an action is an activity,
//! a file an entity and the software an agent; this module does not claim conformance
//! to a profile it has not been validated against.
//!
//! It carries no clock and no host name, so the same audit of the same bytes writes the
//! same document. It says how far an answer moves when the listed decisions change. It
//! does not say which choice is right.

use crate::json::Json;
use crate::{AnalysisAudit, AuditPlan, AuditedRun};
use molframe_core::contract::Provenance;

const CONTEXT: &str = "https://w3id.org/ro/crate/1.1/context";
const SPECIFICATION: &str = "https://w3id.org/ro/crate/1.1";
const DISCLAIMER: &str = "This certificate measures how far the answer moves when the listed \
decisions change. It does not say which choice is correct, and a decision not listed was not \
varied.";

/// The RO-Crate metadata document certifying `audit`, as deterministic JSON.
///
/// Every universe of the plan appears as an action carrying the provenance of the run
/// that answered it. Inputs are listed once, with the SHA-256 of their bytes where the
/// reader recorded it; an input without a digest is listed without one, and the
/// findings count how many there are, since byte-identity cannot be claimed for them.
#[must_use]
pub fn certificate<R: AuditedRun>(plan: &AuditPlan, audit: &AnalysisAudit<R>) -> String {
    let records: Vec<Option<&Provenance>> = audit
        .runs
        .iter()
        .map(|run| run.result.provenance())
        .collect();
    let inputs = Inputs::of(&records);
    let mut graph = vec![descriptor(), root(&inputs)];
    graph.extend(inputs.entities());
    graph.extend(software(&records));
    graph.extend(decisions(plan, &records));
    for (index, run) in audit.runs.iter().enumerate() {
        graph.push(universe(
            index,
            records[index],
            run.result.answer().is_some(),
            &inputs,
        ));
    }
    graph.push(organise(audit, &inputs));
    graph.push(findings(plan, audit, &inputs));
    Json::Object(vec![
        ("@context", Json::text(CONTEXT)),
        ("@graph", Json::Array(graph)),
    ])
    .render()
}

fn object(members: Vec<(&'static str, Json)>) -> Json {
    Json::Object(members)
}

fn property(name: &str, value: Json) -> Json {
    object(vec![
        ("@type", Json::text("PropertyValue")),
        ("name", Json::text(name)),
        ("value", value),
    ])
}

fn descriptor() -> Json {
    object(vec![
        ("@id", Json::text("ro-crate-metadata.json")),
        ("@type", Json::text("CreativeWork")),
        ("conformsTo", Json::reference(SPECIFICATION)),
        ("about", Json::reference("./")),
    ])
}

fn root(inputs: &Inputs) -> Json {
    let mut parts: Vec<Json> = inputs.ids().map(|id| Json::reference(&id)).collect();
    parts.sort_by_key(Json::render);
    object(vec![
        ("@id", Json::text("./")),
        ("@type", Json::text("Dataset")),
        ("name", Json::text("Semantic robustness certificate")),
        ("description", Json::text(DISCLAIMER)),
        ("hasPart", Json::Array(parts)),
        ("mentions", Json::Array(vec![Json::reference("#audit")])),
    ])
}

/// The distinct inputs the runs were computed from.
struct Inputs {
    /// Source, and the SHA-256 when one was recorded.
    entries: Vec<(String, Option<String>)>,
    /// How many runs carried no record of their input at all.
    unrecorded: usize,
}

impl Inputs {
    fn of(records: &[Option<&Provenance>]) -> Self {
        let mut entries: Vec<(String, Option<String>)> = Vec::new();
        let mut unrecorded = 0;
        for record in records {
            let Some(record) = record else {
                unrecorded += 1;
                continue;
            };
            let source = record.input_source.to_string();
            let digest = record.input_digest.map(|digest| digest.to_string());
            if !entries
                .iter()
                .any(|entry| *entry == (source.clone(), digest.clone()))
            {
                entries.push((source, digest));
            }
        }
        Self {
            entries,
            unrecorded,
        }
    }

    fn id(position: usize) -> String {
        format!("#input-{position}")
    }

    fn ids(&self) -> impl Iterator<Item = String> + '_ {
        (0..self.entries.len()).map(Self::id)
    }

    fn without_digest(&self) -> usize {
        self.entries
            .iter()
            .filter(|(_, digest)| digest.is_none())
            .count()
            + self.unrecorded
    }

    fn entities(&self) -> Vec<Json> {
        self.entries
            .iter()
            .enumerate()
            .map(|(position, (source, digest))| {
                let mut members = vec![
                    ("@id", Json::text(Self::id(position))),
                    ("@type", Json::text("File")),
                    ("name", Json::text(source.clone())),
                ];
                if let Some(digest) = digest {
                    members.push((
                        "additionalProperty",
                        property("sha256", Json::text(digest.clone())),
                    ));
                }
                object(members)
            })
            .collect()
    }
}

fn software(records: &[Option<&Provenance>]) -> Vec<Json> {
    let first = records.iter().flatten().next();
    let mut entities = vec![object(vec![
        ("@id", Json::text("#molframe")),
        ("@type", Json::text("SoftwareApplication")),
        ("name", Json::text("molframe")),
        (
            "version",
            first.map_or(Json::Null, |record| Json::text(record.molframe_version)),
        ),
    ])];
    if let Some(algorithm) = first.and_then(|record| record.algorithm.as_ref()) {
        let estimand = first.and_then(|record| record.estimand());
        entities.push(object(vec![
            ("@id", Json::text("#algorithm")),
            ("@type", Json::text("SoftwareApplication")),
            ("name", Json::text(algorithm.name())),
            ("version", Json::text(algorithm.version())),
            ("description", estimand.map_or(Json::Null, Json::text)),
        ]));
    }
    entities
}

fn decisions(plan: &AuditPlan, records: &[Option<&Provenance>]) -> Vec<Json> {
    plan.decisions()
        .iter()
        .map(|decision| {
            let key = format!("policy.{}", decision.field.name());
            let mut values: Vec<String> = Vec::new();
            for record in records.iter().flatten() {
                for (name, value) in record.entries() {
                    if name == key && !values.contains(&value) {
                        values.push(value);
                    }
                }
            }
            object(vec![
                (
                    "@id",
                    Json::text(format!("#decision-{}", decision.field.name())),
                ),
                ("@type", Json::text("PropertyValue")),
                ("name", Json::text(decision.field.name())),
                (
                    "value",
                    Json::Array(values.into_iter().map(Json::Text).collect()),
                ),
                ("description", Json::text(decision.rationale.as_ref())),
                (
                    "additionalProperty",
                    Json::Array(vec![
                        property("uncertainty_class", Json::text(decision.class.name())),
                        property("evidence", Json::text(decision.evidence.as_ref())),
                    ]),
                ),
            ])
        })
        .collect()
}

fn universe(index: usize, record: Option<&Provenance>, answered: bool, inputs: &Inputs) -> Json {
    let mut properties = vec![property(
        "outcome",
        Json::text(if answered {
            "determinate"
        } else {
            "indeterminate"
        }),
    )];
    if let Some(record) = record {
        properties.extend(
            record
                .entries()
                .into_iter()
                .map(|(name, value)| property(&name, Json::Text(value))),
        );
    }
    object(vec![
        ("@id", Json::text(format!("#universe-{index}"))),
        ("@type", Json::text("CreateAction")),
        ("name", Json::text(format!("universe {index}"))),
        ("actionStatus", Json::text("CompletedActionStatus")),
        (
            "instrument",
            Json::Array(vec![
                Json::reference("#molframe"),
                Json::reference("#algorithm"),
            ]),
        ),
        (
            "object",
            Json::Array(inputs.ids().map(|id| Json::reference(&id)).collect()),
        ),
        ("additionalProperty", Json::Array(properties)),
    ])
}

fn organise<R>(audit: &AnalysisAudit<R>, inputs: &Inputs) -> Json {
    let mut objects: Vec<Json> = (0..audit.runs.len())
        .map(|index| Json::reference(&format!("#universe-{index}")))
        .collect();
    objects.extend(inputs.ids().map(|id| Json::reference(&id)));
    object(vec![
        ("@id", Json::text("#audit")),
        ("@type", Json::text("OrganizeAction")),
        ("name", Json::text("semantic robustness audit")),
        ("actionStatus", Json::text("CompletedActionStatus")),
        ("instrument", Json::reference("#molframe")),
        ("object", Json::Array(objects)),
        ("result", Json::reference("#findings")),
    ])
}

fn findings<R>(plan: &AuditPlan, audit: &AnalysisAudit<R>, inputs: &Inputs) -> Json {
    let count = |value: usize| Json::Number(crate::numeric::usize_to_f64(value));
    let mut measured = vec![
        property("metric", Json::text(audit.metric)),
        property("universes", count(audit.runs.len())),
        property("plan_balanced", Json::Bool(plan.balanced())),
        property("plan_skipped_combinations", count(plan.skipped())),
        property(
            "indeterminate_fraction",
            Json::Number(audit.indeterminate_fraction()),
        ),
        property("inputs_without_sha256", count(inputs.without_digest())),
        property(
            "agreement_with_first",
            audit
                .agreement_with_first()
                .map_or(Json::Null, Json::Number),
        ),
    ];
    if let Some(split) = &audit.decomposition {
        measured.push(property(
            "total_variation",
            Json::Number(split.total_variation),
        ));
        measured.push(property("mean_distance", Json::Number(split.mean_distance)));
        measured.push(property("max_distance", Json::Number(split.max_distance)));
        for effect in &split.main_effects {
            let name = effect.field.name();
            measured.push(property(
                &format!("main_effect.{name}.mean_change"),
                Json::Number(effect.mean_change),
            ));
            measured.push(property(
                &format!("main_effect.{name}.share"),
                Json::Number(effect.share),
            ));
        }
        for pair in &split.interactions {
            measured.push(property(
                &format!("interaction.{}.{}", pair.first.name(), pair.second.name()),
                Json::Number(pair.share),
            ));
        }
        measured.push(property("higher_order", Json::Number(split.higher_order)));
        for share in &split.shapley {
            measured.push(property(
                &format!("shapley.{}", share.key.name()),
                Json::Number(share.share),
            ));
        }
        for share in &split.by_class {
            measured.push(property(
                &format!("class.{}", share.key.name()),
                Json::Number(share.share),
            ));
        }
    }
    object(vec![
        ("@id", Json::text("#findings")),
        ("@type", Json::text("Dataset")),
        ("name", Json::text("measured robustness")),
        ("description", Json::text(DISCLAIMER)),
        ("variableMeasured", Json::Array(measured)),
    ])
}

#[cfg(test)]
#[path = "certificate_tests.rs"]
mod tests;
