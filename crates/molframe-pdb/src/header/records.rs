use molframe_core::Structure;

/// Stable extension key for ordered legacy PDB metadata records.
pub const PDB_HEADERS_EXTENSION: &str = "molframe.pdb.headers.v1";

/// One fixed-column metadata record exactly as deposited.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PdbHeaderRecord {
    name: Box<str>,
    line: Box<str>,
}

impl PdbHeaderRecord {
    pub(crate) fn new(name: &str, line: &str) -> Self {
        Self {
            name: name.into(),
            line: line.into(),
        }
    }

    /// Six-column record name, trimmed.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Original fixed-column line without its line ending.
    #[must_use]
    pub fn line(&self) -> &str {
        &self.line
    }
}

/// Every metadata record in archive order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PdbHeaders {
    records: Vec<PdbHeaderRecord>,
}

impl PdbHeaders {
    pub(crate) fn push(&mut self, name: &str, line: &str) {
        self.records.push(PdbHeaderRecord::new(name, line));
    }

    /// Whether no metadata records were deposited.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// All metadata records in archive order.
    #[must_use]
    pub fn records(&self) -> &[PdbHeaderRecord] {
        &self.records
    }

    /// Records of one type, retaining their archive order.
    pub fn named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a PdbHeaderRecord> {
        self.records
            .iter()
            .filter(move |record| record.name.as_ref() == name)
    }

    /// Classification from columns 11–50 of `HEADER`.
    #[must_use]
    pub fn classification(&self) -> Option<&str> {
        self.header_field(10, 50)
    }

    /// Deposition date from columns 51–59 of `HEADER`.
    #[must_use]
    pub fn deposition_date(&self) -> Option<&str> {
        self.header_field(50, 59)
    }

    fn header_field(&self, start: usize, end: usize) -> Option<&str> {
        let header = self.named("HEADER").next()?.line();
        let value = header.get(start..end)?.trim();
        (!value.is_empty()).then_some(value)
    }
}

/// Access to preserved legacy metadata.
pub trait PdbHeadersExt {
    /// Metadata attached by the legacy PDB reader.
    fn pdb_headers(&self) -> Option<&PdbHeaders>;
}

impl PdbHeadersExt for Structure {
    fn pdb_headers(&self) -> Option<&PdbHeaders> {
        self.extensions().get(PDB_HEADERS_EXTENSION)
    }
}

/// Whether a record belongs to a metadata section of PDB format 3.3.
pub(crate) fn is_metadata_record(name: &str) -> bool {
    matches!(
        name,
        "HEADER"
            | "OBSLTE"
            | "TITLE"
            | "SPLIT"
            | "CAVEAT"
            | "COMPND"
            | "SOURCE"
            | "KEYWDS"
            | "EXPDTA"
            | "NUMMDL"
            | "MDLTYP"
            | "AUTHOR"
            | "REVDAT"
            | "SPRSDE"
            | "JRNL"
            | "REMARK"
            | "DBREF"
            | "DBREF1"
            | "DBREF2"
            | "SEQADV"
            | "SEQRES"
            | "MODRES"
            | "HET"
            | "HETNAM"
            | "HETSYN"
            | "FORMUL"
            | "HELIX"
            | "SHEET"
            | "SSBOND"
            | "LINK"
            | "CISPEP"
            | "SITE"
            | "CRYST1"
            | "ORIGX1"
            | "ORIGX2"
            | "ORIGX3"
            | "SCALE1"
            | "SCALE2"
            | "SCALE3"
            | "MTRIX1"
            | "MTRIX2"
            | "MTRIX3"
    )
}
