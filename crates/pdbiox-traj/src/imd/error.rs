//! IMD connection and protocol errors.

/// A transport, framing, validation, or protocol-negotiation failure.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ImdError {
    /// The network or supplied stream failed.
    #[error("IMD transport error: {0}")]
    Io(#[from] std::io::Error),
    /// The peer did not begin with a supported handshake.
    #[error("invalid or unsupported IMD handshake")]
    InvalidHandshake,
    /// A header carried an unknown message type.
    #[error("unknown IMD message type {0}")]
    UnknownMessage(i32),
    /// A message length contradicted its type or configured resource limits.
    #[error("invalid IMD message length {length} for type {message_type}")]
    InvalidLength {
        /// Numeric wire message type.
        message_type: i32,
        /// Header length field.
        length: i32,
    },
    /// A coordinate, energy, force, rate, or atom index was invalid.
    #[error("IMD message contains an invalid numeric value")]
    InvalidValue,
    /// An operation was attempted after a terminal peer message.
    #[error("IMD session is closed")]
    Closed,
}
