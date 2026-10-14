//! Final assembly of parsed frames, topology and format extensions.

use super::ReadState;
use super::secondary;
use crate::header::PDB_HEADERS_EXTENSION;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::ReadResult;
use molframe_core::structure::{CoordinateStore, Structure};

impl ReadState<'_> {
    pub(in crate::reader) fn finish_anisou(&mut self) {
        if self.anisou.is_empty() {
            return;
        }
        let mut builder = molframe_core::AnisotropyTableBuilder::new();
        for (serial, u) in self.anisou.drain(..) {
            let Some(atom) = self.serial_to_atom.get(&serial).copied() else {
                continue;
            };
            builder.push(molframe_core::AnisotropicDisplacement { atom, u });
        }
        self.data.anisotropy = builder.finish();
    }

    pub(in crate::reader) fn finish(mut self) -> ReadResult {
        self.close_chain();
        self.verify_frame_len();
        if self.later_models.values {
            self.findings
                .push(Diagnostic::new(Code::W3015).with_message(
                    "occupancies or B-factors differ between models; model 1's values are kept",
                ));
        }
        if !self.saw_atoms {
            // A row that looks like a coordinate record but whose field layout
            // does not parse is a different failure from a file that holds no
            // coordinates, and saying which is the difference between a fixable
            // report and a baffling one.
            let finding = if self.malformed_atom > 0 {
                Diagnostic::new(Code::E1201)
                    .with_message("coordinate records are present but their fields do not parse")
                    .with_context("rows", self.malformed_atom.to_string())
            } else {
                Diagnostic::new(Code::E1001).with_message("the file contains no coordinate records")
            };
            self.findings.push(finding);
            return Err(self.findings.finish());
        }

        self.finish_bonds();
        self.finish_anisou();
        let (states, sources) = secondary::read(&self.data, &self.headers);
        self.data.secondary_structure = states.into();
        self.data.secondary_source = sources.into();
        self.finish_variant_annotations();
        if self.topology_locked {
            self.finish_dense_frames();
        } else {
            let builder = std::mem::take(&mut self.builder);
            let (chunks, coords) = builder.finish();
            self.data.chunks = chunks.into();
            self.frames.push(coords);
        }
        if let Some(resolution) = self
            .headers
            .resolution()
            .and_then(|value| num_traits::ToPrimitive::to_f32(&value))
        {
            self.data.entry.resolution = Some(resolution);
        }
        if !self.headers.is_empty() {
            self.data
                .extensions
                .insert(PDB_HEADERS_EXTENSION, self.headers);
        }
        let mut chains = 0..0;
        if let Ok(count) = u32::try_from(self.data.topology.chains.len()) {
            chains = 0..count;
        } else {
            self.findings.push(
                Diagnostic::new(Code::E1901)
                    .with_message("the chain table exceeds the supported index range"),
            );
        }
        for number in self.model_numbers {
            if self
                .data
                .topology
                .models
                .push(number, chains.clone())
                .is_err()
            {
                self.findings.push(Diagnostic::new(Code::E3001));
            }
        }
        self.data.coords = match self.frames.len() {
            1 => match self.frames.pop() {
                Some(block) => CoordinateStore::Single(block),
                None => CoordinateStore::Single(CoordinateBlock::new()),
            },
            _ => CoordinateStore::Dense {
                frames: self.frames,
            },
        };
        let structure = Structure::new(self.data);
        self.options.finish(structure, self.findings.finish())
    }
}
