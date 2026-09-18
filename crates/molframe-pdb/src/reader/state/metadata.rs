use super::ReadState;
use crate::fixed;
use crate::header::is_metadata_record;
use crate::reader::lines::Line;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::span::ByteSpan;
use molframe_core::structure::UnitCell;

impl ReadState<'_> {
    pub(super) fn observe_metadata(&mut self, record: &str, line: &Line<'_>) {
        if is_metadata_record(record) {
            self.headers.push(record, line.text);
        }
    }

    pub(super) fn cell(&mut self, line: &Line<'_>) {
        let lengths = [
            fixed::real(line.text, 7, 15),
            fixed::real(line.text, 16, 24),
            fixed::real(line.text, 25, 33),
        ];
        let angles = [
            fixed::real(line.text, 34, 40),
            fixed::real(line.text, 41, 47),
            fixed::real(line.text, 48, 54),
        ];
        let (Some(a), Some(b), Some(c)) = (lengths[0], lengths[1], lengths[2]) else {
            self.findings.push(
                Diagnostic::new(Code::E1202)
                    .with_message("cell lengths could not be read")
                    .at(ByteSpan::empty(line.at)),
            );
            return;
        };
        let (Some(alpha), Some(beta), Some(gamma)) = (angles[0], angles[1], angles[2]) else {
            return;
        };
        self.data.cell = Some(UnitCell {
            lengths: [a, b, c],
            angles: [alpha, beta, gamma],
        });
    }

    pub(super) fn header(&mut self, line: &Line<'_>) {
        let id = fixed::text(line.text, 63, 66);
        if !id.is_empty() {
            self.data.entry.id = Some(id.into());
        }
    }

    pub(super) fn title(&mut self, line: &Line<'_>) {
        append(&mut self.data.entry.title, fixed::text(line.text, 11, 80));
    }

    pub(super) fn method(&mut self, line: &Line<'_>) {
        append(&mut self.data.entry.method, fixed::text(line.text, 11, 79));
    }
}

fn append(target: &mut Option<Box<str>>, text: &str) {
    if text.is_empty() {
        return;
    }
    *target = Some(match target.take() {
        Some(existing) => format!("{existing} {text}").into(),
        None => text.into(),
    });
}
