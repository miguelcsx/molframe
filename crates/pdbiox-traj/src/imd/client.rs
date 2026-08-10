//! Stateful IMD client over TCP or a caller-supplied duplex stream.

use std::io::{Read, Write};
use std::net::{TcpStream, ToSocketAddrs};

use super::ImdError;
use super::model::{ImdConnectionOptions, ImdForce, ImdLimits, ImdMessage, ImdPeerEndian};
use super::wire::{
    DISCONNECT, GO, KILL, PAUSE, TRATE, negotiate, read_message, write_control, write_forces,
};

/// A negotiated IMD v2 session.
///
/// Reads are exact and streaming: at most one complete message is resident in
/// memory. The configured atom limit is checked before any peer-sized
/// allocation.
#[derive(Debug)]
pub struct ImdClient<S> {
    stream: S,
    peer_endian: ImdPeerEndian,
    limits: ImdLimits,
    closed: bool,
}

impl ImdClient<TcpStream> {
    /// Connects to a simulation with default socket and allocation policy.
    ///
    /// # Errors
    ///
    /// Returns a transport error or rejects an invalid IMD handshake.
    pub fn connect(address: impl ToSocketAddrs) -> Result<Self, ImdError> {
        Self::connect_with_options(address, ImdConnectionOptions::default())
    }

    /// Connects with explicit socket and allocation policy.
    ///
    /// # Errors
    ///
    /// Returns a transport error, invalid option, or handshake failure.
    pub fn connect_with_options(
        address: impl ToSocketAddrs,
        options: ImdConnectionOptions,
    ) -> Result<Self, ImdError> {
        validate_limits(options.limits)?;
        let stream = TcpStream::connect(address)?;
        stream.set_read_timeout(options.read_timeout)?;
        stream.set_write_timeout(options.write_timeout)?;
        stream.set_nodelay(options.no_delay)?;
        Self::from_stream(stream, options.limits)
    }
}

impl<S: Read + Write> ImdClient<S> {
    /// Negotiates IMD v2 over an already-connected duplex stream.
    ///
    /// This is useful for Unix sockets, TLS wrappers, test transports, and
    /// embedding. The simulation must speak first with its handshake; the
    /// client replies with `GO` only after version and byte order are known.
    ///
    /// # Errors
    ///
    /// Returns a transport error, invalid limit, or handshake failure.
    pub fn from_stream(mut stream: S, limits: ImdLimits) -> Result<Self, ImdError> {
        validate_limits(limits)?;
        let peer_endian = negotiate(&mut stream)?;
        write_control(&mut stream, GO, 0)?;
        Ok(Self {
            stream,
            peer_endian,
            limits,
            closed: false,
        })
    }

    /// Receives exactly one complete peer message.
    ///
    /// # Errors
    ///
    /// Returns a transport or protocol error, including premature EOF.
    pub fn receive(&mut self) -> Result<ImdMessage, ImdError> {
        self.ensure_open()?;
        let message = read_message(&mut self.stream, self.peer_endian, self.limits)?;
        if matches!(message, ImdMessage::Disconnect | ImdMessage::Kill) {
            self.closed = true;
        }
        Ok(message)
    }

    /// Sends atom-targeted steering forces in canonical units.
    ///
    /// The wire conversion to kcal mol⁻¹ Å⁻¹ and the peer's byte order is
    /// performed once while writing the contiguous protocol blocks.
    ///
    /// # Errors
    ///
    /// Returns an error for a closed session, invalid force/index, allocation
    /// policy violation, or transport failure.
    pub fn send_forces(&mut self, forces: &[ImdForce]) -> Result<(), ImdError> {
        self.ensure_open()?;
        write_forces(&mut self.stream, self.peer_endian, forces, self.limits)
    }

    /// Sends the protocol's pause-toggle command.
    ///
    /// # Errors
    ///
    /// Returns an error when the session is closed or writing fails.
    pub fn toggle_pause(&mut self) -> Result<(), ImdError> {
        self.send_control(PAUSE, 0)
    }

    /// Sends `GO`, which starts or resumes a simulation.
    ///
    /// # Errors
    ///
    /// Returns an error when the session is closed or writing fails.
    pub fn resume(&mut self) -> Result<(), ImdError> {
        self.send_control(GO, 0)
    }

    /// Requests a positive coordinate transmission interval.
    ///
    /// # Errors
    ///
    /// Returns an error for zero, values outside signed 32-bit range, a closed
    /// session, or transport failure.
    pub fn set_transmission_rate(&mut self, steps: u32) -> Result<(), ImdError> {
        let rate = i32::try_from(steps).map_err(|_| ImdError::InvalidValue)?;
        if rate == 0 {
            return Err(ImdError::InvalidValue);
        }
        self.send_control(TRATE, rate)
    }

    /// Detaches while leaving the simulation running.
    ///
    /// # Errors
    ///
    /// Returns an error when the session is already closed or writing fails.
    pub fn disconnect(&mut self) -> Result<(), ImdError> {
        self.send_control(DISCONNECT, 0)?;
        self.closed = true;
        Ok(())
    }

    /// Requests termination of the simulation job.
    ///
    /// # Errors
    ///
    /// Returns an error when the session is already closed or writing fails.
    pub fn kill(&mut self) -> Result<(), ImdError> {
        self.send_control(KILL, 0)?;
        self.closed = true;
        Ok(())
    }

    /// Returns the negotiated peer payload byte order.
    #[must_use]
    pub const fn peer_endian(&self) -> ImdPeerEndian {
        self.peer_endian
    }

    /// Whether a terminal local or peer command closed the session.
    #[must_use]
    pub const fn is_closed(&self) -> bool {
        self.closed
    }

    /// Consumes the client and returns its transport.
    #[must_use]
    pub fn into_inner(self) -> S {
        self.stream
    }

    fn send_control(&mut self, message_type: i32, length: i32) -> Result<(), ImdError> {
        self.ensure_open()?;
        write_control(&mut self.stream, message_type, length)
    }

    const fn ensure_open(&self) -> Result<(), ImdError> {
        if self.closed {
            Err(ImdError::Closed)
        } else {
            Ok(())
        }
    }
}

fn validate_limits(limits: ImdLimits) -> Result<(), ImdError> {
    if limits.max_atoms == 0 || i32::try_from(limits.max_atoms).is_err() {
        Err(ImdError::InvalidValue)
    } else {
        Ok(())
    }
}
