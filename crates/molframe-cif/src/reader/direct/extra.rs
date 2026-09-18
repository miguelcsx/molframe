//! Multiplexing an external event projection into the direct structure pass.

use super::stream::{DirectOutput, DirectSink};
use crate::document::CifValueRef;
use crate::lexer::Quoting;
use crate::parser::{CifEventSink, CifScalar, ColumnId, ValueSink};
use molframe_core::io::ReadOptions;
use molframe_core::span::ByteSpan;

pub(super) struct ProjectedSink<'options, 'input, S> {
    direct: DirectSink<'options, 'input>,
    projection: S,
}

impl<'options, S> ProjectedSink<'options, '_, S> {
    pub(super) fn new(
        options: &'options ReadOptions,
        keep_category: fn(&str) -> bool,
        input_bytes: usize,
        projection: S,
    ) -> Self {
        Self {
            direct: DirectSink::new(options, keep_category, input_bytes),
            projection,
        }
    }
}

impl<'input, S> ValueSink<'input> for ProjectedSink<'_, 'input, S>
where
    S: CifEventSink,
{
    type Output = (DirectOutput, S::Output);

    fn block(&mut self, name: &str) {
        self.direct.block(name);
        self.projection.block(name);
    }

    fn begin_loop(&mut self, tags: &[(Box<str>, Box<str>)]) {
        self.direct.begin_loop(tags);
        self.projection.begin_loop(tags);
    }

    fn value(
        &mut self,
        column: ColumnId,
        category: &str,
        item: &str,
        text: &'input str,
        quoting: Quoting,
        span: ByteSpan,
    ) {
        self.direct
            .value(column, category, item, text, quoting, span);
        if self.projection.accepts_category(category) {
            self.projection.value(
                category,
                item,
                CifScalar::from(CifValueRef::parse(text, quoting)),
                span,
            );
        }
    }

    fn end_row(&mut self) {
        self.direct.end_row();
        self.projection.end_row();
    }

    fn end_loop(&mut self) {
        self.direct.end_loop();
        self.projection.end_loop();
    }

    fn finish(self) -> Self::Output {
        (self.direct.finish(), self.projection.finish())
    }
}
