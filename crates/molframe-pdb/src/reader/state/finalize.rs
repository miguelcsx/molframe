//! Final assembly of parsed frames, topology and format extensions.

use super::ReadState;
use crate::header::PDB_HEADERS_EXTENSION;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::ReadResult;
use molframe_core::structure::{CoordinateStore, Structure};

impl ReadState<'_> {
    pub(in crate::reader) fn finish(mut self) -> ReadResult {
        self.close_chain();
        self.verify_frame_len();
        if !self.saw_atoms {
            self.findings.push(
                Diagnostic::new(Code::E1001)
                    .with_message("the file contains no coordinate records"),
            );
            return Err(self.findings.finish());
        }

        self.finish_bonds();
        self.finish_variant_annotations();
        if self.topology_locked {
            self.finish_dense_frames();
        } else {
            let builder = std::mem::take(&mut self.builder);
            let (chunks, coords) = builder.finish();
            self.data.chunks = chunks.into();
            self.frames.push(coords);
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
