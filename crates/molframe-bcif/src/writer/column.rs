//! Compact typed column builders for the structure projection.

use crate::container::EncodedColumn;
use crate::encode::StringColumnEncoder;
use crate::{encode_floats, encode_integers};
use molframe_cif::CanonicalValue;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::element::Element;

pub(super) struct TextColumnBuilder {
    name: &'static str,
    data: StringColumnEncoder,
    mask: MaskBuilder,
}

impl TextColumnBuilder {
    pub(super) fn new(name: &'static str, rows: usize) -> Self {
        Self {
            name,
            data: StringColumnEncoder::with_capacity(rows),
            mask: MaskBuilder::new(rows),
        }
    }

    pub(super) fn push(&mut self, value: CanonicalValue<&str>) {
        match value {
            CanonicalValue::Present(text) => {
                self.data.push(text);
                self.mask.push(0);
            }
            CanonicalValue::Inapplicable => {
                self.data.push("");
                self.mask.push(1);
            }
            CanonicalValue::Unknown => {
                self.data.push("");
                self.mask.push(2);
            }
        }
    }

    pub(super) fn push_element(&mut self, value: CanonicalValue<Element>) {
        match value {
            CanonicalValue::Present(element) => {
                self.data.push_ascii_uppercase(element.symbol());
                self.mask.push(0);
            }
            CanonicalValue::Inapplicable => {
                self.data.push("");
                self.mask.push(1);
            }
            CanonicalValue::Unknown => {
                self.data.push("");
                self.mask.push(2);
            }
        }
    }

    pub(super) fn finish(self) -> Result<EncodedColumn, Diagnostic> {
        let Self { name, data, mask } = self;
        let data = data.finish()?;
        Ok(EncodedColumn {
            name: name.to_owned(),
            data,
            mask: mask.finish()?,
        })
    }
}

pub(super) struct IntegerColumnBuilder {
    name: &'static str,
    data: Vec<i64>,
    mask: MaskBuilder,
}

impl IntegerColumnBuilder {
    pub(super) fn new(name: &'static str, rows: usize) -> Self {
        Self {
            name,
            data: Vec::with_capacity(rows),
            mask: MaskBuilder::new(rows),
        }
    }

    pub(super) fn push(&mut self, value: CanonicalValue<i64>) {
        match value {
            CanonicalValue::Present(number) => {
                self.data.push(number);
                self.mask.push(0);
            }
            CanonicalValue::Inapplicable => {
                self.data.push(0);
                self.mask.push(1);
            }
            CanonicalValue::Unknown => {
                self.data.push(0);
                self.mask.push(2);
            }
        }
    }

    pub(super) fn finish(self) -> Result<EncodedColumn, Diagnostic> {
        let Self { name, data, mask } = self;
        let encoded = encode_integers(&data)?;
        drop(data);
        Ok(EncodedColumn {
            name: name.to_owned(),
            data: encoded,
            mask: mask.finish()?,
        })
    }
}

pub(super) struct FloatColumnBuilder {
    name: &'static str,
    data: Vec<f64>,
    mask: MaskBuilder,
}

impl FloatColumnBuilder {
    pub(super) fn new(name: &'static str, rows: usize) -> Self {
        Self {
            name,
            data: Vec::with_capacity(rows),
            mask: MaskBuilder::new(rows),
        }
    }

    pub(super) fn push(&mut self, value: CanonicalValue<f64>) {
        match value {
            CanonicalValue::Present(number) => {
                self.data.push(number);
                self.mask.push(0);
            }
            CanonicalValue::Inapplicable => {
                self.data.push(0.0);
                self.mask.push(1);
            }
            CanonicalValue::Unknown => {
                self.data.push(0.0);
                self.mask.push(2);
            }
        }
    }

    pub(super) fn finish(self) -> Result<EncodedColumn, Diagnostic> {
        let Self { name, data, mask } = self;
        let encoded = encode_floats(&data)?;
        drop(data);
        Ok(EncodedColumn {
            name: name.to_owned(),
            data: encoded,
            mask: mask.finish()?,
        })
    }
}

struct MaskBuilder {
    rows: usize,
    length: usize,
    values: Option<Vec<i64>>,
}

impl MaskBuilder {
    const fn new(rows: usize) -> Self {
        Self {
            rows,
            length: 0,
            values: None,
        }
    }

    fn push(&mut self, value: i64) {
        if value == 0 {
            if let Some(values) = &mut self.values {
                values.push(0);
            }
        } else if let Some(values) = &mut self.values {
            values.push(value);
        } else {
            let mut values = Vec::with_capacity(self.rows);
            values.resize(self.length, 0);
            values.push(value);
            self.values = Some(values);
        }
        self.length += 1;
    }

    fn finish(self) -> Result<Option<crate::EncodedData>, Diagnostic> {
        match self.values {
            Some(values) => encode_integers(&values).map(Some),
            None => Ok(None),
        }
    }
}
