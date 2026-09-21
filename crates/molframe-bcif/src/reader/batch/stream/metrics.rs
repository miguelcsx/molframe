//! Retained dictionary metrics used for exact pre-consumption reservations.

use super::{ColumnStream, TextStream, ValueStream};

impl ColumnStream {
    pub(in crate::reader::batch) fn max_text_width(&self) -> usize {
        match &self.values {
            ValueStream::Text(stream) => stream.max_width(),
            ValueStream::Integer(_) | ValueStream::Float(_) => 0,
        }
    }
}

impl TextStream {
    fn max_width(&self) -> usize {
        let maximum = self
            .offsets
            .windows(2)
            .filter_map(|pair| pair[1].checked_sub(pair[0]))
            .filter_map(|width| usize::try_from(width).ok())
            .max();
        let Some(width) = maximum else {
            return 0;
        };
        width
    }
}
