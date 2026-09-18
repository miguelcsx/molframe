//! Allocation-free typed views over compact quality-assessment categories.

use crate::{CompactColumn, ModelCategory, ModelCif};

/// Definition shared by quality metric value categories.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MetricDefinition<'a> {
    /// Dictionary metric identifier.
    pub id: &'a str,
    /// Human-readable name.
    pub name: Option<&'a str>,
    /// Controlled metric type, such as `pLDDT`, `PAE`, or `pTM`.
    pub metric_type: &'a str,
    /// Controlled mode: global, local, or local-pairwise.
    pub mode: &'a str,
    /// Software group responsible for the value.
    pub software_group_id: Option<&'a str>,
}

/// One model-level metric value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GlobalMetric<'a> {
    /// `ModelCIF` model identifier.
    pub model_id: &'a str,
    /// Definition identifier.
    pub metric_id: &'a str,
    /// Metric value.
    pub value: f64,
}

/// One residue-level metric value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalMetric<'a> {
    /// `ModelCIF` model identifier.
    pub model_id: &'a str,
    /// Label chain identifier.
    pub chain_id: &'a str,
    /// Label residue number.
    pub sequence_id: i64,
    /// Component identifier when supplied.
    pub component_id: Option<&'a str>,
    /// Definition identifier.
    pub metric_id: &'a str,
    /// Metric value.
    pub value: f64,
}

/// One residue-pair metric value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PairwiseMetric<'a> {
    /// `ModelCIF` model identifier.
    pub model_id: &'a str,
    /// First label chain identifier.
    pub first_chain_id: &'a str,
    /// First label residue number.
    pub first_sequence_id: i64,
    /// Second label chain identifier.
    pub second_chain_id: &'a str,
    /// Second label residue number.
    pub second_sequence_id: i64,
    /// Definition identifier.
    pub metric_id: &'a str,
    /// Metric value.
    pub value: f64,
}

/// Typed views over all compact `ModelCIF` QA definitions and values.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QualityMetrics<'a> {
    model: &'a ModelCif,
}

impl<'a> QualityMetrics<'a> {
    pub(crate) const fn new(model: &'a ModelCif) -> Self {
        Self { model }
    }

    /// Metric definitions in source-row order.
    #[must_use]
    pub fn definitions(self) -> MetricDefinitions<'a> {
        MetricDefinitions::new(self.model.category("ma_qa_metric"))
    }

    /// Model-level metrics in source-row order.
    #[must_use]
    pub fn global(self) -> GlobalMetrics<'a> {
        GlobalMetrics::new(self.model.category("ma_qa_metric_global"))
    }

    /// Residue-level metrics in source-row order.
    #[must_use]
    pub fn local(self) -> LocalMetrics<'a> {
        LocalMetrics::new(self.model.category("ma_qa_metric_local"))
    }

    /// Residue-pair metrics in source-row order.
    #[must_use]
    pub fn pairwise(self) -> PairwiseMetrics<'a> {
        PairwiseMetrics::new(self.model.category("ma_qa_metric_local_pairwise"))
    }

    /// Residue-level pLDDT values, including dictionary variants whose type begins with pLDDT.
    pub fn plddt(self) -> impl Iterator<Item = LocalMetric<'a>> + 'a {
        self.local().filter(move |value| {
            self.definition(value.metric_id)
                .is_some_and(|definition| definition.metric_type.starts_with("pLDDT"))
        })
    }

    /// Residue-pair PAE values.
    pub fn pae(self) -> impl Iterator<Item = PairwiseMetric<'a>> + 'a {
        self.pairwise().filter(move |value| {
            self.definition(value.metric_id)
                .is_some_and(|definition| definition.metric_type.eq_ignore_ascii_case("PAE"))
        })
    }

    /// Model-level pTM values.
    pub fn ptm(self) -> impl Iterator<Item = GlobalMetric<'a>> + 'a {
        self.global().filter(move |value| {
            self.definition(value.metric_id)
                .is_some_and(|definition| definition.metric_type.eq_ignore_ascii_case("pTM"))
        })
    }

    fn definition(self, id: &str) -> Option<MetricDefinition<'a>> {
        self.definitions().find(|definition| definition.id == id)
    }
}

#[derive(Clone, Copy, Debug)]
struct DefinitionColumns<'a> {
    id: &'a CompactColumn,
    name: Option<&'a CompactColumn>,
    metric_type: &'a CompactColumn,
    mode: &'a CompactColumn,
    software_group_id: Option<&'a CompactColumn>,
}

/// Iterator over metric-definition rows.
#[derive(Clone, Debug)]
pub struct MetricDefinitions<'a> {
    columns: Option<DefinitionColumns<'a>>,
    row: usize,
    rows: usize,
}

impl<'a> MetricDefinitions<'a> {
    fn new(category: Option<&'a ModelCategory>) -> Self {
        let columns = category.and_then(|category| {
            Some(DefinitionColumns {
                id: category.column("id")?,
                name: category.column("name"),
                metric_type: category.column("type")?,
                mode: category.column("mode")?,
                software_group_id: category.column("software_group_id"),
            })
        });
        Self {
            columns,
            row: 0,
            rows: category.map_or(0, ModelCategory::row_count),
        }
    }
}

impl<'a> Iterator for MetricDefinitions<'a> {
    type Item = MetricDefinition<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let columns = self.columns?;
        while self.row < self.rows {
            let row = self.row;
            self.row += 1;
            if let Some(parsed) = definition_at(columns, row) {
                return Some(parsed);
            }
        }
        None
    }
}

#[derive(Clone, Copy, Debug)]
struct GlobalColumns<'a> {
    model_id: &'a CompactColumn,
    metric_id: &'a CompactColumn,
    value: &'a CompactColumn,
}

/// Iterator over model-level metric rows.
#[derive(Clone, Debug)]
pub struct GlobalMetrics<'a> {
    columns: Option<GlobalColumns<'a>>,
    row: usize,
    rows: usize,
}

impl<'a> GlobalMetrics<'a> {
    fn new(category: Option<&'a ModelCategory>) -> Self {
        let columns = category.and_then(|category| {
            Some(GlobalColumns {
                model_id: category.column("model_id")?,
                metric_id: category.column("metric_id")?,
                value: category.column("metric_value")?,
            })
        });
        Self {
            columns,
            row: 0,
            rows: category.map_or(0, ModelCategory::row_count),
        }
    }
}

impl<'a> Iterator for GlobalMetrics<'a> {
    type Item = GlobalMetric<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let columns = self.columns?;
        while self.row < self.rows {
            let row = self.row;
            self.row += 1;
            if let Some(parsed) = global_at(columns, row) {
                return Some(parsed);
            }
        }
        None
    }
}

#[derive(Clone, Copy, Debug)]
struct LocalColumns<'a> {
    model_id: &'a CompactColumn,
    chain_id: &'a CompactColumn,
    sequence_id: &'a CompactColumn,
    component_id: Option<&'a CompactColumn>,
    metric_id: &'a CompactColumn,
    value: &'a CompactColumn,
}

/// Iterator over residue-level metric rows.
#[derive(Clone, Debug)]
pub struct LocalMetrics<'a> {
    columns: Option<LocalColumns<'a>>,
    row: usize,
    rows: usize,
}

impl<'a> LocalMetrics<'a> {
    fn new(category: Option<&'a ModelCategory>) -> Self {
        let columns = category.and_then(|category| {
            Some(LocalColumns {
                model_id: category.column("model_id")?,
                chain_id: category.column("label_asym_id")?,
                sequence_id: category.column("label_seq_id")?,
                component_id: category.column("label_comp_id"),
                metric_id: category.column("metric_id")?,
                value: category.column("metric_value")?,
            })
        });
        Self {
            columns,
            row: 0,
            rows: category.map_or(0, ModelCategory::row_count),
        }
    }
}

impl<'a> Iterator for LocalMetrics<'a> {
    type Item = LocalMetric<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let columns = self.columns?;
        while self.row < self.rows {
            let row = self.row;
            self.row += 1;
            if let Some(parsed) = local_at(columns, row) {
                return Some(parsed);
            }
        }
        None
    }
}

#[derive(Clone, Copy, Debug)]
struct PairwiseColumns<'a> {
    model_id: &'a CompactColumn,
    first_chain_id: &'a CompactColumn,
    first_sequence_id: &'a CompactColumn,
    second_chain_id: &'a CompactColumn,
    second_sequence_id: &'a CompactColumn,
    metric_id: &'a CompactColumn,
    value: &'a CompactColumn,
}

/// Iterator over residue-pair metric rows.
#[derive(Clone, Debug)]
pub struct PairwiseMetrics<'a> {
    columns: Option<PairwiseColumns<'a>>,
    row: usize,
    rows: usize,
}

impl<'a> PairwiseMetrics<'a> {
    fn new(category: Option<&'a ModelCategory>) -> Self {
        let columns = category.and_then(|category| {
            Some(PairwiseColumns {
                model_id: category.column("model_id")?,
                first_chain_id: category.column("label_asym_id_1")?,
                first_sequence_id: category.column("label_seq_id_1")?,
                second_chain_id: category.column("label_asym_id_2")?,
                second_sequence_id: category.column("label_seq_id_2")?,
                metric_id: category.column("metric_id")?,
                value: category.column("metric_value")?,
            })
        });
        Self {
            columns,
            row: 0,
            rows: category.map_or(0, ModelCategory::row_count),
        }
    }
}

impl<'a> Iterator for PairwiseMetrics<'a> {
    type Item = PairwiseMetric<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let columns = self.columns?;
        while self.row < self.rows {
            let row = self.row;
            self.row += 1;
            if let Some(parsed) = pairwise_at(columns, row) {
                return Some(parsed);
            }
        }
        None
    }
}

fn optional_id(column: Option<&CompactColumn>, row: usize) -> Option<&str> {
    column.and_then(|column| column.identifier(row))
}

fn definition_at(columns: DefinitionColumns<'_>, row: usize) -> Option<MetricDefinition<'_>> {
    Some(MetricDefinition {
        id: columns.id.identifier(row)?,
        name: optional_id(columns.name, row),
        metric_type: columns.metric_type.identifier(row)?,
        mode: columns.mode.identifier(row)?,
        software_group_id: optional_id(columns.software_group_id, row),
    })
}

fn global_at(columns: GlobalColumns<'_>, row: usize) -> Option<GlobalMetric<'_>> {
    let value = columns.value.value(row)?.number()?;
    if !value.is_finite() {
        return None;
    }
    Some(GlobalMetric {
        model_id: columns.model_id.identifier(row)?,
        metric_id: columns.metric_id.identifier(row)?,
        value,
    })
}

fn local_at(columns: LocalColumns<'_>, row: usize) -> Option<LocalMetric<'_>> {
    let value = columns.value.value(row)?.number()?;
    if !value.is_finite() {
        return None;
    }
    Some(LocalMetric {
        model_id: columns.model_id.identifier(row)?,
        chain_id: columns.chain_id.identifier(row)?,
        sequence_id: columns.sequence_id.value(row)?.integer()?,
        component_id: optional_id(columns.component_id, row),
        metric_id: columns.metric_id.identifier(row)?,
        value,
    })
}

fn pairwise_at(columns: PairwiseColumns<'_>, row: usize) -> Option<PairwiseMetric<'_>> {
    let value = columns.value.value(row)?.number()?;
    if !value.is_finite() {
        return None;
    }
    Some(PairwiseMetric {
        model_id: columns.model_id.identifier(row)?,
        first_chain_id: columns.first_chain_id.identifier(row)?,
        first_sequence_id: columns.first_sequence_id.value(row)?.integer()?,
        second_chain_id: columns.second_chain_id.identifier(row)?,
        second_sequence_id: columns.second_sequence_id.value(row)?.integer()?,
        metric_id: columns.metric_id.identifier(row)?,
        value,
    })
}
