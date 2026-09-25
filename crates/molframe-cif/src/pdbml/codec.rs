//! PDBML/XML conversion through the lossless mmCIF document model.

use crate::document::{CifValue, DataBlock, Document};
use crate::lexer::Quoting;
use molframe_core::ByteSpan;
use molframe_core::io::ReadOptions;
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};
use quick_xml::{Reader, Writer};

/// Malformed or unsupported PDBML framing.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PdbmlError {
    /// XML tokenization or entity decoding failed.
    #[error("invalid PDBML/XML: {0}")]
    Xml(#[from] quick_xml::Error),
    /// XML output could not be written.
    #[error("could not write PDBML/XML: {0}")]
    Io(#[from] std::io::Error),
    /// The document has no PDBML datablock.
    #[error("PDBML/XML has no datablock")]
    MissingDatablock,
    /// PDBML rows/categories are nested inconsistently.
    #[error("invalid PDBML category or row nesting")]
    InvalidNesting,
    /// XML output was not UTF-8, which cannot represent PDBML.
    #[error("PDBML writer produced invalid UTF-8")]
    InvalidUtf8,
    /// An XML entity reference was not predefined or numeric.
    #[error("unsupported PDBML XML entity reference")]
    InvalidEntity,
}

/// Failure either in XML conversion or existing mmCIF semantic lowering.
#[derive(Debug)]
#[non_exhaustive]
pub enum PdbmlReadError {
    /// PDBML syntax/framing failure.
    Pdbml(PdbmlError),
    /// Normal structure-lowering findings that prevented a structure.
    Findings(Vec<molframe_core::Diagnostic>),
}

impl std::fmt::Display for PdbmlReadError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Pdbml(error) => error.fmt(formatter),
            Self::Findings(findings) => {
                write!(
                    formatter,
                    "PDBML semantic lowering failed with {} findings",
                    findings.len()
                )
            }
        }
    }
}
impl std::error::Error for PdbmlReadError {}

/// Parses PDBML/XML into the same lossless document used by mmCIF.
/// Unknown categories and items are retained in source order.
///
/// # Errors
///
/// Returns an error for malformed XML, absent datablocks, or inconsistent
/// category/row nesting.
pub fn parse_pdbml_document(bytes: &[u8]) -> Result<Document, PdbmlError> {
    let mut reader = Reader::from_reader(bytes);
    reader.config_mut().trim_text(false);
    let mut document = Document::new();
    let mut block = None;
    let mut category: Option<String> = None;
    let mut row: Option<Row> = None;
    let mut item: Option<Item> = None;
    loop {
        match reader.read_event()? {
            Event::Start(element) => {
                let name = local_name(element.name().as_ref()).to_owned();
                if name == "datablock" {
                    block = Some(DataBlock::new(datablock_name(&element, &reader)?));
                } else if let Some(value) = name.strip_suffix("Category") {
                    if category.is_some() || row.is_some() {
                        return Err(PdbmlError::InvalidNesting);
                    }
                    category = Some(value.to_owned());
                } else if category.as_deref() == Some(&name) {
                    if row.is_some() {
                        return Err(PdbmlError::InvalidNesting);
                    }
                    row = Some(Row::from_attributes(&element, &reader)?);
                } else if row.is_some() {
                    item = Some(Item {
                        name,
                        text: String::new(),
                        nil: nil_attribute(&element, &reader)?,
                    });
                }
            }
            Event::Empty(element) if row.is_some() => {
                let name = local_name(element.name().as_ref()).to_owned();
                let value = if nil_attribute(&element, &reader)? {
                    CifValue::Unknown
                } else {
                    CifValue::Text("".into())
                };
                if let Some(active) = &mut row {
                    active.values.push((name, value));
                }
            }
            Event::Text(text) => {
                if let Some(active) = &mut item {
                    active
                        .text
                        .push_str(&text.decode().map_err(|_| PdbmlError::InvalidUtf8)?);
                }
            }
            Event::CData(text) => {
                if let Some(active) = &mut item {
                    active
                        .text
                        .push_str(&text.decode().map_err(|_| PdbmlError::InvalidUtf8)?);
                }
            }
            Event::GeneralRef(reference) => {
                if let Some(active) = &mut item {
                    let reference = reference.decode().map_err(|_| PdbmlError::InvalidUtf8)?;
                    active.text.push_str(&resolve_reference(&reference)?);
                }
            }
            Event::End(element) => {
                let name = local_name(element.name().as_ref()).to_owned();
                if item.as_ref().is_some_and(|active| active.name == name) {
                    let active = item.take().ok_or(PdbmlError::InvalidNesting)?;
                    let value = if active.nil {
                        CifValue::Unknown
                    } else {
                        CifValue::parse(active.text.trim(), Quoting::Bare)
                    };
                    row.as_mut()
                        .ok_or(PdbmlError::InvalidNesting)?
                        .values
                        .push((active.name, value));
                } else if category.as_deref() == Some(&name) {
                    let active = row.take().ok_or(PdbmlError::InvalidNesting)?;
                    let block = block.as_mut().ok_or(PdbmlError::MissingDatablock)?;
                    append_row(block, &name, active);
                } else if name.ends_with("Category") {
                    category = None;
                } else if name == "datablock" {
                    document.push(block.take().ok_or(PdbmlError::MissingDatablock)?);
                }
            }
            Event::Eof => break,
            _ => {}
        }
    }
    if document.is_empty() {
        return Err(PdbmlError::MissingDatablock);
    }
    Ok(document)
}

/// Parses and semantically lowers PDBML through the existing mmCIF lowering.
///
/// # Errors
///
/// Returns the XML error or the full semantic finding list.
pub fn read_pdbml(
    bytes: &[u8],
    options: &ReadOptions,
) -> Result<
    (
        Document,
        molframe_core::Structure,
        Vec<molframe_core::Diagnostic>,
    ),
    PdbmlReadError,
> {
    let document = parse_pdbml_document(bytes).map_err(PdbmlReadError::Pdbml)?;
    let (structure, findings) =
        crate::lower(&document, options).map_err(PdbmlReadError::Findings)?;
    options
        .finish(structure, findings)
        .map(|(structure, findings)| (document, structure, findings))
        .map_err(PdbmlReadError::Findings)
}

/// Writes one lossless document as PDBML/XML.
///
/// Both CIF missing sentinels map to XML Schema `nil`, because PDBML has one
/// missing representation rather than CIF's two distinct sentinels.
///
/// # Errors
///
/// Returns an error when the document has other than one datablock or XML
/// generation fails.
pub fn write_pdbml(document: &Document) -> Result<String, PdbmlError> {
    if document.len() != 1 {
        return Err(PdbmlError::MissingDatablock);
    }
    let block = document.first_block().ok_or(PdbmlError::MissingDatablock)?;
    let mut writer = Writer::new_with_indent(Vec::new(), b' ', 3);
    writer.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None)))?;
    let mut root = BytesStart::new("PDBx:datablock");
    root.push_attribute(("datablockName", block.name()));
    root.push_attribute(("xmlns:PDBx", "http://pdbml.pdb.org/schema/pdbx-v50.xsd"));
    root.push_attribute(("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance"));
    writer.write_event(Event::Start(root))?;
    for category in block.categories() {
        let category_name = format!("PDBx:{}Category", category.name());
        writer.write_event(Event::Start(BytesStart::new(&category_name)))?;
        for row in 0..category.row_count() {
            let row_name = format!("PDBx:{}", category.name());
            writer.write_event(Event::Start(BytesStart::new(&row_name)))?;
            for item_name in category.items() {
                let value = match category.value(item_name, row) {
                    Some(value) => value,
                    None => &CifValue::Unknown,
                };
                write_item(&mut writer, item_name, value)?;
            }
            writer.write_event(Event::End(BytesEnd::new(&row_name)))?;
        }
        writer.write_event(Event::End(BytesEnd::new(&category_name)))?;
    }
    writer.write_event(Event::End(BytesEnd::new("PDBx:datablock")))?;
    String::from_utf8(writer.into_inner()).map_err(|_| PdbmlError::InvalidUtf8)
}

struct Row {
    values: Vec<(String, CifValue)>,
}
impl Row {
    fn from_attributes(
        element: &BytesStart<'_>,
        reader: &Reader<&[u8]>,
    ) -> Result<Self, PdbmlError> {
        let mut values = Vec::new();
        for attribute in element.attributes().with_checks(false) {
            let attribute = attribute.map_err(quick_xml::Error::from)?;
            let name = local_name(attribute.key.as_ref());
            if name == "nil" || name == "schemaLocation" || name.starts_with("xmlns") {
                continue;
            }
            let text = attribute.decoded_and_normalized_value(
                quick_xml::XmlVersion::Implicit1_0,
                reader.decoder(),
            )?;
            values.push((name.to_owned(), CifValue::parse(&text, Quoting::Bare)));
        }
        Ok(Self { values })
    }
}
struct Item {
    name: String,
    text: String,
    nil: bool,
}

fn datablock_name(element: &BytesStart<'_>, reader: &Reader<&[u8]>) -> Result<String, PdbmlError> {
    let _ = reader;

    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(quick_xml::Error::from)?;

        if local_name(attribute.key.as_ref()) == "datablockName" {
            return Ok(attribute
                .normalized_value(quick_xml::XmlVersion::Implicit1_0)?
                .into_owned());
        }
    }

    Err(PdbmlError::MissingDatablock)
}

fn nil_attribute(element: &BytesStart<'_>, reader: &Reader<&[u8]>) -> Result<bool, PdbmlError> {
    let _ = reader;

    for attribute in element.attributes().with_checks(false) {
        let attribute = attribute.map_err(quick_xml::Error::from)?;

        if local_name(attribute.key.as_ref()) == "nil" {
            return Ok(attribute
                .normalized_value(quick_xml::XmlVersion::Implicit1_0)?
                .eq_ignore_ascii_case("true"));
        }
    }

    Ok(false)
}

fn local_name(name: &[u8]) -> &str {
    let local = name
        .iter()
        .rposition(|byte| *byte == b':')
        .map_or(name, |position| &name[position + 1..]);
    std::str::from_utf8(local).map_or("", |value| value)
}
fn resolve_reference(reference: &str) -> Result<String, PdbmlError> {
    if let Some(value) = quick_xml::escape::resolve_predefined_entity(reference) {
        return Ok(value.to_owned());
    }
    let number = if let Some(hexadecimal) = reference.strip_prefix("#x") {
        u32::from_str_radix(hexadecimal, 16).ok()
    } else {
        reference
            .strip_prefix('#')
            .and_then(|value| value.parse().ok())
    };
    number
        .and_then(char::from_u32)
        .map(|value| value.to_string())
        .ok_or(PdbmlError::InvalidEntity)
}
fn append_row(block: &mut DataBlock, category_name: &str, row: Row) {
    let category = block.category_mut(category_name, ByteSpan::default());
    let row_index = category.row_count();
    let existing: Vec<_> = category.items().map(str::to_owned).collect();
    for item in existing {
        let value = row
            .values
            .iter()
            .find(|(name, _)| name == &item)
            .map_or(CifValue::Unknown, |(_, value)| value.clone());
        category.column_mut(&item).push(value, Quoting::Bare);
    }
    for (item, value) in row.values {
        if category.column(&item).is_some() {
            continue;
        }
        let column = category.column_mut(&item);
        for _ in 0..row_index {
            column.push(CifValue::Unknown, Quoting::Bare);
        }
        column.push(value, Quoting::Bare);
    }
}
fn write_item(
    writer: &mut Writer<Vec<u8>>,
    item_name: &str,
    value: &CifValue,
) -> Result<(), PdbmlError> {
    let name = format!("PDBx:{item_name}");
    match value {
        CifValue::Unknown | CifValue::Inapplicable => {
            let mut element = BytesStart::new(&name);
            element.push_attribute(("xsi:nil", "true"));
            writer.write_event(Event::Empty(element))?;
        }
        CifValue::Text(text) => write_text_element(writer, &name, text)?,
        CifValue::Integer(value) => write_text_element(writer, &name, &value.to_string())?,
        CifValue::Float(value) => write_text_element(writer, &name, &value.to_string())?,
    }
    Ok(())
}
fn write_text_element(
    writer: &mut Writer<Vec<u8>>,
    name: &str,
    text: &str,
) -> Result<(), PdbmlError> {
    writer.write_event(Event::Start(BytesStart::new(name)))?;
    writer.write_event(Event::Text(BytesText::new(text)))?;
    writer.write_event(Event::End(BytesEnd::new(name)))?;
    Ok(())
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
