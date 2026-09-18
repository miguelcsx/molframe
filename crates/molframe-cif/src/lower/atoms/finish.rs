//! Final assembly of lowered coordinate frames and model topology.

use super::AtomBuilder;
use molframe_core::coords::CoordinateBlock;
use molframe_core::diagnostic::{Code, Diagnostic, Diagnostics};
use molframe_core::structure::{CoordinateStore, StructureData};

impl AtomBuilder<'_> {
    pub(crate) fn finish(mut self) -> (StructureData, Diagnostics, CoordinateStore) {
        self.close_residue();
        self.close_chain();
        self.verify_frame_len();
        let (chunks, coords) = self.builder.finish();
        if self.data.chunks.is_empty() {
            self.data.chunks = chunks.into();
        }
        self.frames.push(coords);

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

        let coords = match self.frames.len() {
            1 => match self.frames.pop() {
                Some(block) => CoordinateStore::Single(block),
                None => CoordinateStore::Single(CoordinateBlock::new()),
            },
            _ => CoordinateStore::Dense {
                frames: self.frames,
            },
        };
        (self.data, self.findings, coords)
    }

    pub(crate) fn abort(self) -> Diagnostics {
        self.findings
    }
}
