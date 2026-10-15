//! Lazily lowered CCD provider.
//!
//! Scans the dictionary bytes once for `data_` headers (linear in input size)
//! and lowers only the blocks a caller asks for, caching each result.

use super::ComponentProvider;
use super::lowering::component;
use crate::model::Component;
use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::InputBuffer;
use std::sync::Arc;

/// A component dictionary that lowers a block only when it is asked for.
///
/// `read_ccd` parses and lowers every block of a dictionary even when a
/// structure uses a few dozen components. This provider scans the raw bytes
/// once for `data_` headers, then parses and lowers just the requested blocks
/// and caches them. A malformed block is reported when it is first requested,
/// not at open time, so use `read_ccd` when whole-dictionary validation is the
/// goal.
pub struct LazyCifProvider {
    input: InputBuffer,
    version: DictionaryVersion,
    index: std::collections::HashMap<Box<str>, std::ops::Range<usize>>,
    cache: std::sync::RwLock<std::collections::HashMap<Box<str>, Option<Arc<Component>>>>,
}

impl std::fmt::Debug for LazyCifProvider {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("LazyCifProvider")
            .field("version", &self.version)
            .field("blocks", &self.index.len())
            .field("lowered", &self.lowered_count())
            .finish_non_exhaustive()
    }
}

impl LazyCifProvider {
    /// Indexes the component blocks of a dictionary without lowering them.
    #[must_use]
    pub fn new(input: InputBuffer, version: DictionaryVersion) -> Self {
        let index = index_blocks(input.as_bytes());
        Self {
            input,
            version,
            index,
            cache: std::sync::RwLock::new(std::collections::HashMap::new()),
        }
    }

    /// Number of indexed component blocks.
    #[must_use]
    pub fn block_count(&self) -> usize {
        self.index.len()
    }

    /// Number of blocks lowered so far.
    #[must_use]
    pub fn lowered_count(&self) -> usize {
        match self.cache.read() {
            Ok(cache) => cache.len(),
            Err(_) => 0,
        }
    }

    fn lower_block(
        &self,
        id: &str,
        range: &std::ops::Range<usize>,
    ) -> Result<Component, Diagnostic> {
        let slice = match self.input.as_bytes().get(range.clone()) {
            Some(slice) => slice,
            None => &[],
        };
        let buffer = InputBuffer::from_bytes(slice.to_vec());
        let (document, _) = molframe_cif::parse(&buffer).map_err(first_finding)?;
        let block = document.blocks().next().ok_or_else(|| {
            Diagnostic::new(Code::E2001)
                .with_context("block", id)
                .in_category("chem_comp")
        })?;
        let component = component(block).map_err(first_finding)?;
        if component.id.as_ref() == id {
            Ok(component)
        } else {
            Err(Diagnostic::new(Code::E2005)
                .in_category("chem_comp")
                .with_context("id", component.id.to_string())
                .with_context("block", id))
        }
    }
}

fn first_finding(findings: Vec<Diagnostic>) -> Diagnostic {
    match findings.into_iter().next() {
        Some(finding) => finding,
        None => Diagnostic::new(Code::E2001).in_category("chem_comp"),
    }
}

/// Byte ranges of every `data_` block, keyed by block name.
fn index_blocks(bytes: &[u8]) -> std::collections::HashMap<Box<str>, std::ops::Range<usize>> {
    let mut starts: Vec<(usize, Box<str>)> = Vec::new();
    let mut line = 0;
    while line < bytes.len() {
        let rest = match bytes.get(line..) {
            Some(rest) => rest,
            None => &[],
        };
        let end = rest.iter().position(|byte| *byte == b'\n');
        if rest.starts_with(b"data_") {
            let header_end = match end {
                Some(end) => end,
                None => rest.len(),
            };
            let header = &rest[..header_end];
            let name = &header[5..];
            let name = match name.split(u8::is_ascii_whitespace).next() {
                Some(name) => name,
                None => &[],
            };
            if let Ok(name) = std::str::from_utf8(name) {
                starts.push((line, name.into()));
            }
        }
        let Some(offset) = end else {
            break;
        };
        line += offset + 1;
    }
    let mut index = std::collections::HashMap::with_capacity(starts.len());
    for (position, (start, name)) in starts.iter().enumerate() {
        let stop = match starts.get(position + 1) {
            Some(next) => next.0,
            None => bytes.len(),
        };
        index.entry(name.clone()).or_insert(*start..stop);
    }
    index
}

impl ComponentProvider for LazyCifProvider {
    fn get(&self, component_id: &str) -> Result<Option<Arc<Component>>, Diagnostic> {
        let cached = self
            .cache
            .read()
            .ok()
            .and_then(|cache| cache.get(component_id).cloned());
        if let Some(hit) = cached {
            return Ok(hit);
        }
        let lowered = match self.index.get(component_id) {
            Some(range) => Some(Arc::new(self.lower_block(component_id, range)?)),
            None => None,
        };
        if let Ok(mut cache) = self.cache.write() {
            cache.insert(component_id.into(), lowered.clone());
        }
        Ok(lowered)
    }

    fn version(&self) -> &DictionaryVersion {
        &self.version
    }
}
