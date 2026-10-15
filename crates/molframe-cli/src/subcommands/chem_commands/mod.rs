//! Small-molecule chemistry commands over the dictionary-backed kernels.

mod dispatch;
mod peoe;
mod smarts;

pub(crate) use dispatch::execute;
