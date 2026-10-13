//! `BinaryCIF` write verbs.

use molframe_core::diagnostic::Findings;

use crate::structure::Structure;

/// Renders deterministic `BinaryCIF` bytes in memory.
///
/// # Errors
///
/// Returns a diagnostic if a projected column cannot be represented.
pub fn write_bcif(structure: &Structure) -> Result<Vec<u8>, Findings> {
    write_bcif_with_options(structure, &molframe_cif::CifWriteOptions::new())
}

/// Renders deterministic `BinaryCIF` in memory with explicit identifier decisions.
///
/// # Errors
///
/// Returns a diagnostic if canonical preflight or binary encoding fails.
pub fn write_bcif_with_options(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<Vec<u8>, Findings> {
    write_bcif_with_findings(structure, options).map(|(bytes, _)| bytes)
}

/// Renders `BinaryCIF` like [`write_bcif_with_options`] and also reports what
/// the write rewrote or could not carry.
///
/// A renamed chain or shifted residue number is `W4102`; an extension with no
/// `BinaryCIF` counterpart is `W4101`. Lossless input returns no findings.
///
/// # Errors
///
/// Returns a diagnostic if canonical preflight or binary encoding fails.
pub fn write_bcif_with_findings(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<(Vec<u8>, Findings), Findings> {
    let bytes = encode(structure, options)?;
    Ok((bytes, super::cif_family_findings(structure)))
}

fn encode(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<Vec<u8>, Findings> {
    let options = super::complete_cif_options(structure.engine(), options);
    #[cfg(feature = "modelcif")]
    if let Some(model) = structure
        .engine()
        .extensions()
        .get::<molframe_modelcif::ModelCif>(molframe_modelcif::MODEL_CIF_EXTENSION)
    {
        return write_with_model(structure, model, &options);
    }
    molframe_bcif::write_structure_with_options(structure.engine(), &options)
        .map_err(Findings::from)
}

/// Appends the retained `ModelCIF` tables, encoded column by column.
///
/// Cells are borrowed from the compact columns and the dense PAE matrix; no
/// text and no parsed document is built.
#[cfg(feature = "modelcif")]
fn write_with_model(
    structure: &Structure,
    model: &molframe_modelcif::ModelCif,
    options: &molframe_cif::CifWriteOptions,
) -> Result<Vec<u8>, Findings> {
    struct Table<'a>(molframe_modelcif::ExportTable<'a>);

    impl molframe_bcif::ColumnTable for Table<'_> {
        fn name(&self) -> &str {
            self.0.name()
        }
        fn row_count(&self) -> usize {
            self.0.row_count()
        }
        fn column_count(&self) -> usize {
            self.0.column_count()
        }
        fn item(&self, column: usize) -> &str {
            self.0.item(column)
        }
        fn cell(&self, column: usize, row: usize) -> molframe_bcif::Cell<'_> {
            use molframe_bcif::Cell;
            use molframe_modelcif::ExportCell;
            match self.0.cell(column, row) {
                ExportCell::Inapplicable => Cell::Inapplicable,
                ExportCell::Unknown => Cell::Unknown,
                ExportCell::Text(value) => Cell::Text(value),
                ExportCell::Integer(value) => Cell::Integer(value),
                ExportCell::Float(value) => Cell::Float(value),
            }
        }
    }

    let tables: Vec<Table<'_>> = model.export_tables().into_iter().map(Table).collect();
    let views: Vec<&dyn molframe_bcif::ColumnTable> = tables
        .iter()
        .map(|table| table as &dyn molframe_bcif::ColumnTable)
        .collect();
    molframe_bcif::write_structure_with_tables(structure.engine(), options, &views)
        .map_err(Findings::from)
}
