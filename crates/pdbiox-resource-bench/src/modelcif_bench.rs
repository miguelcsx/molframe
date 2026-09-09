//! Process-isolated compact `ModelCIF` projection for an external official corpus.

use super::{ResourceRecord, measure_retained_case};
use std::hint::black_box;
use std::path::Path;

pub(super) fn read_file(path: &Path) -> Result<ResourceRecord, String> {
    let bytes = std::fs::read(path)
        .map_err(|error| format!("{} cannot be read: {error}", path.display()))?;
    let input = pdbiox::InputBuffer::from_bytes(bytes);
    measure_retained_case("modelcif_file", || {
        let (model, findings) = pdbiox::modelcif::read_compact(&input)
            .map_err(|error| format!("{} projection failed: {error}", path.display()))?;
        if !findings.is_empty() {
            return Err(format!(
                "{} projection reported findings: {findings:?}",
                path.display()
            ));
        }
        let rows = model
            .category("ma_qa_metric_local_pairwise")
            .map_or(0, pdbiox::ModelCategory::row_count);
        let digest = u64::try_from(rows).map_err(|_| "ModelCIF rows exceed u64".to_owned())?;
        black_box(&model);
        Ok((digest, None, model))
    })
}
