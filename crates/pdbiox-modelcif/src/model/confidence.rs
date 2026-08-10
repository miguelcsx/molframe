/// Definition shared by quality metric value categories.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetricDefinition {
    /// Dictionary metric identifier.
    pub id: Box<str>,
    /// Human-readable name.
    pub name: Option<Box<str>>,
    /// Controlled metric type, such as `pLDDT`, `PAE`, or `pTM`.
    pub metric_type: Box<str>,
    /// Controlled mode: global, local, or local-pairwise.
    pub mode: Box<str>,
    /// Software group responsible for the value.
    pub software_group_id: Option<Box<str>>,
}

/// One model-level metric value.
#[derive(Clone, Debug, PartialEq)]
pub struct GlobalMetric {
    /// `ModelCIF` model identifier.
    pub model_id: Box<str>,
    /// Definition identifier.
    pub metric_id: Box<str>,
    /// Metric value.
    pub value: f64,
}

/// One residue-level metric value.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalMetric {
    /// `ModelCIF` model identifier.
    pub model_id: Box<str>,
    /// Label chain identifier.
    pub chain_id: Box<str>,
    /// Label residue number.
    pub sequence_id: i64,
    /// Component identifier when supplied.
    pub component_id: Option<Box<str>>,
    /// Definition identifier.
    pub metric_id: Box<str>,
    /// Metric value.
    pub value: f64,
}

/// One residue-pair metric value.
#[derive(Clone, Debug, PartialEq)]
pub struct PairwiseMetric {
    /// `ModelCIF` model identifier.
    pub model_id: Box<str>,
    /// First label chain identifier.
    pub first_chain_id: Box<str>,
    /// First label residue number.
    pub first_sequence_id: i64,
    /// Second label chain identifier.
    pub second_chain_id: Box<str>,
    /// Second label residue number.
    pub second_sequence_id: i64,
    /// Definition identifier.
    pub metric_id: Box<str>,
    /// Metric value.
    pub value: f64,
}

/// All typed `ModelCIF` QA definitions and values.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct QualityMetrics {
    /// Metric definitions.
    pub definitions: Vec<MetricDefinition>,
    /// Model-level values.
    pub global: Vec<GlobalMetric>,
    /// Residue-level values.
    pub local: Vec<LocalMetric>,
    /// Residue-pair values.
    pub pairwise: Vec<PairwiseMetric>,
}

impl QualityMetrics {
    /// Residue-level pLDDT values, including dictionary variants whose type begins with pLDDT.
    pub fn plddt(&self) -> impl Iterator<Item = &LocalMetric> {
        self.local.iter().filter(|value| {
            self.definition(&value.metric_id)
                .is_some_and(|definition| definition.metric_type.starts_with("pLDDT"))
        })
    }

    /// Residue-pair PAE values.
    pub fn pae(&self) -> impl Iterator<Item = &PairwiseMetric> {
        self.pairwise.iter().filter(|value| {
            self.definition(&value.metric_id)
                .is_some_and(|definition| definition.metric_type.eq_ignore_ascii_case("PAE"))
        })
    }

    /// Model-level pTM values.
    pub fn ptm(&self) -> impl Iterator<Item = &GlobalMetric> {
        self.global.iter().filter(|value| {
            self.definition(&value.metric_id)
                .is_some_and(|definition| definition.metric_type.eq_ignore_ascii_case("pTM"))
        })
    }

    fn definition(&self, id: &str) -> Option<&MetricDefinition> {
        self.definitions
            .iter()
            .find(|definition| definition.id.as_ref() == id)
    }
}
