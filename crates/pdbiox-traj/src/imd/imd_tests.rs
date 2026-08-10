use std::io::{Read, Write};
use std::net::TcpListener;
use std::thread;

use super::{
    ImdClient, ImdConnectionOptions, ImdError, ImdForce, ImdLimits, ImdMessage, ImdPeerEndian,
};
use crate::imd::wire::{
    DISCONNECT, ENERGIES, FCOORDS, GO, HANDSHAKE, HEADER_BYTES, IMD_VERSION, MDCOMM, PAUSE, TRATE,
};

const KCAL_TO_KJ: f32 = 4.184;

struct FragmentedStream {
    input: Vec<u8>,
    offset: usize,
    output: Vec<u8>,
    fragment: usize,
}

impl FragmentedStream {
    fn new(input: Vec<u8>, fragment: usize) -> Self {
        Self {
            input,
            offset: 0,
            output: Vec::new(),
            fragment,
        }
    }
}

impl Read for FragmentedStream {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if self.offset == self.input.len() {
            return Ok(0);
        }
        let available = self.input.len() - self.offset;
        let count = available.min(buffer.len()).min(self.fragment);
        buffer[..count].copy_from_slice(&self.input[self.offset..self.offset + count]);
        self.offset += count;
        Ok(count)
    }
}

impl Write for FragmentedStream {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.output.extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

#[test]
fn fragmented_stream_negotiates_and_receives_coordinates_and_energies() {
    let endian = opposite_native_endian();
    let mut input = handshake(endian);
    input.extend(header(FCOORDS, 2));
    for value in [1.0_f32, 2.0, 3.0, 4.0, 5.0, 6.0] {
        input.extend(float_bytes(value, endian));
    }
    input.extend(header(ENERGIES, 1));
    input.extend(int_bytes(42, endian));
    for value in [300.0_f32, 10.0, 8.0, 1.0, 2.0, 3.0, 0.5, 0.25, 0.125] {
        input.extend(float_bytes(value, endian));
    }
    let stream = FragmentedStream::new(input, 1);
    let mut client = ImdClient::from_stream(stream, ImdLimits { max_atoms: 2 })
        .expect("negotiate fragmented stream");
    assert_eq!(client.peer_endian(), endian);
    assert_eq!(
        client.receive().expect("coordinates"),
        ImdMessage::Coordinates(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0],])
    );
    let ImdMessage::Energies(energies) = client.receive().expect("energies") else {
        panic!("expected energy message");
    };
    assert_eq!(energies.step, 42);
    assert_close(energies.temperature, 300.0);
    assert_close(energies.total, 10.0 * KCAL_TO_KJ);
    assert_close(energies.improper, 0.125 * KCAL_TO_KJ);
    assert_eq!(&client.into_inner().output[..HEADER_BYTES], &header(GO, 0));
}

#[test]
fn steering_forces_use_peer_endian_and_protocol_units() {
    let endian = ImdPeerEndian::Big;
    let stream = FragmentedStream::new(handshake(endian), 8);
    let mut client =
        ImdClient::from_stream(stream, ImdLimits { max_atoms: 4 }).expect("negotiate stream");
    client
        .send_forces(&[
            ImdForce {
                atom: 7,
                force: [4.184, 8.368, 12.552],
            },
            ImdForce {
                atom: 9,
                force: [-4.184, 0.0, 2.092],
            },
        ])
        .expect("send forces");
    let output = client.into_inner().output;
    let payload = &output[HEADER_BYTES..];
    assert_eq!(&payload[..HEADER_BYTES], &header(MDCOMM, 2));
    assert_eq!(&payload[8..12], &int_bytes(7, endian));
    assert_eq!(&payload[12..16], &int_bytes(9, endian));
    assert_close(read_float(&payload[16..20], endian), 1.0);
    assert_close(read_float(&payload[36..40], endian), 0.5);
}

#[test]
fn configured_atom_limit_rejects_header_before_payload_read() {
    let endian = native_endian();
    let mut input = handshake(endian);
    input.extend(header(FCOORDS, 3));
    let stream = FragmentedStream::new(input, 2);
    let mut client =
        ImdClient::from_stream(stream, ImdLimits { max_atoms: 2 }).expect("negotiate stream");
    assert!(matches!(
        client.receive(),
        Err(ImdError::InvalidLength {
            message_type: FCOORDS,
            length: 3
        })
    ));
}

#[test]
fn received_steering_forces_are_converted_and_atom_indexed() {
    let endian = opposite_native_endian();
    let mut input = handshake(endian);
    input.extend(header(MDCOMM, 2));
    input.extend(int_bytes(3, endian));
    input.extend(int_bytes(8, endian));
    for value in [1.0_f32, 2.0, 3.0, -1.0, -2.0, -3.0] {
        input.extend(float_bytes(value, endian));
    }
    let stream = FragmentedStream::new(input, 3);
    let mut client =
        ImdClient::from_stream(stream, ImdLimits { max_atoms: 2 }).expect("negotiate stream");
    let ImdMessage::Forces(forces) = client.receive().expect("force message") else {
        panic!("expected force message");
    };
    assert_eq!(
        forces.iter().map(|force| force.atom).collect::<Vec<_>>(),
        [3, 8]
    );
    assert_close(forces[0].force[2], 3.0 * KCAL_TO_KJ);
    assert_close(forces[1].force[0], -KCAL_TO_KJ);
}

#[test]
fn control_commands_are_framed_and_terminal_state_is_enforced() {
    let stream = FragmentedStream::new(handshake(native_endian()), 8);
    let mut client =
        ImdClient::from_stream(stream, ImdLimits { max_atoms: 1 }).expect("negotiate stream");
    client.toggle_pause().expect("pause");
    client.resume().expect("resume");
    client.set_transmission_rate(25).expect("rate");
    client.disconnect().expect("disconnect");
    assert!(client.is_closed());
    assert!(matches!(client.resume(), Err(ImdError::Closed)));
    let output = client.into_inner().output;
    let expected = [
        header(GO, 0),
        header(PAUSE, 0),
        header(GO, 0),
        header(TRATE, 25),
        header(DISCONNECT, 0),
    ]
    .into_iter()
    .flatten()
    .collect::<Vec<_>>();
    assert_eq!(output, expected);
}

#[test]
fn handshake_rejects_unknown_version_before_sending_go() {
    let mut input = HANDSHAKE.to_be_bytes().to_vec();
    input.extend(int_bytes(IMD_VERSION + 1, native_endian()));
    let stream = FragmentedStream::new(input, 8);
    assert!(matches!(
        ImdClient::from_stream(stream, ImdLimits { max_atoms: 1 }),
        Err(ImdError::InvalidHandshake)
    ));
}

#[test]
fn tcp_connect_performs_handshake_and_go_before_streaming() {
    let listener = TcpListener::bind("127.0.0.1:0").expect("bind loopback");
    let address = listener.local_addr().expect("listener address");
    let server = thread::spawn(move || {
        let (mut stream, _) = listener.accept().expect("accept client");
        stream
            .write_all(&handshake(native_endian()))
            .expect("handshake");
        let mut go = [0u8; HEADER_BYTES];
        stream.read_exact(&mut go).expect("read GO");
        assert_eq!(go, header(GO, 0));
        stream
            .write_all(&header(FCOORDS, 1))
            .expect("coordinate header");
        for value in [7.0_f32, 8.0, 9.0] {
            stream.write_all(&value.to_ne_bytes()).expect("coordinate");
        }
    });
    let mut client = ImdClient::connect_with_options(
        address,
        ImdConnectionOptions {
            limits: ImdLimits { max_atoms: 1 },
            ..ImdConnectionOptions::default()
        },
    )
    .expect("connect client");
    assert_eq!(
        client.receive().expect("coordinate message"),
        ImdMessage::Coordinates(vec![[7.0, 8.0, 9.0]])
    );
    server.join().expect("server thread");
}

fn handshake(endian: ImdPeerEndian) -> Vec<u8> {
    let mut bytes = HANDSHAKE.to_be_bytes().to_vec();
    bytes.extend(int_bytes(IMD_VERSION, endian));
    bytes
}

fn header(message_type: i32, length: i32) -> [u8; HEADER_BYTES] {
    let mut bytes = [0u8; HEADER_BYTES];
    bytes[..4].copy_from_slice(&message_type.to_be_bytes());
    bytes[4..].copy_from_slice(&length.to_be_bytes());
    bytes
}

fn int_bytes(value: i32, endian: ImdPeerEndian) -> [u8; 4] {
    match endian {
        ImdPeerEndian::Little => value.to_le_bytes(),
        ImdPeerEndian::Big => value.to_be_bytes(),
    }
}

fn float_bytes(value: f32, endian: ImdPeerEndian) -> [u8; 4] {
    match endian {
        ImdPeerEndian::Little => value.to_le_bytes(),
        ImdPeerEndian::Big => value.to_be_bytes(),
    }
}

fn read_float(bytes: &[u8], endian: ImdPeerEndian) -> f32 {
    let array = [bytes[0], bytes[1], bytes[2], bytes[3]];
    match endian {
        ImdPeerEndian::Little => f32::from_le_bytes(array),
        ImdPeerEndian::Big => f32::from_be_bytes(array),
    }
}

const fn native_endian() -> ImdPeerEndian {
    if cfg!(target_endian = "little") {
        ImdPeerEndian::Little
    } else {
        ImdPeerEndian::Big
    }
}

const fn opposite_native_endian() -> ImdPeerEndian {
    if cfg!(target_endian = "little") {
        ImdPeerEndian::Big
    } else {
        ImdPeerEndian::Little
    }
}

fn assert_close(actual: f32, expected: f32) {
    assert!((actual - expected).abs() < 1.0e-5);
}
