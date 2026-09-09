//! An mmCIF byte stream generated as it is consumed.
//!
//! This is what makes a hundred-gigabyte measurement possible on a workstation.
//! The stream is produced one batch of coordinate rows at a time from a tile and
//! a seed, so it occupies one buffer no matter how long it runs, touches no
//! disk, and is reproducible from the seed without being stored.
//!
//! The emitted block carries the fifteen `atom_site` items a structural reader
//! actually consumes, which puts a row at about ninety bytes — the same figure
//! the direct reader uses to estimate row counts from input size.
//!
//! Chain labels carry the copy ordinal, so every copy is its own chain. Without
//! that, consecutive copies would present identical residue annotations and the
//! reader would fall back to inferring boundaries from file order — a real path,
//! but a rare one, and measuring it instead of the ordinary path would make the
//! fixture report the wrong thing. It also lets the identifier dictionary grow
//! with structure size, which a fixed set of labels would leave untested.

use super::format::{push_fixed, push_i32, push_u64};
use super::seed::Seed;
use super::tile::{Placement, Tile};
use std::io::{self, Read};

/// Coordinate rows formatted per buffer refill.
///
/// Large enough that the per-refill bookkeeping disappears against the
/// formatting, small enough that the buffer stays inside second-level cache.
const ROWS_PER_REFILL: usize = 512;

/// The block header, emitted once before any coordinate row.
const HEADER: &str = "data_synthetic\n\
loop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_entity_id\n\
_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n\
_atom_site.pdbx_PDB_model_num\n";

/// An mmCIF stream generated on demand from a tile and a seed.
///
/// # Examples
///
/// ```
/// use pdbiox_bench::{Sample, Seed, SyntheticCifSource, Tile, structure};
/// use std::io::Read;
///
/// let tile = Tile::from_structure(&structure(Sample::Tiny), Seed::new(1));
/// let mut source = SyntheticCifSource::new(tile, 4);
/// let mut head = [0_u8; 14];
/// source.read_exact(&mut head)?;
/// assert_eq!(&head, b"data_synthetic");
/// # Ok::<(), std::io::Error>(())
/// ```
#[derive(Debug)]
pub struct SyntheticCifSource {
    tile: Tile,
    copies: u64,
    copy: u64,
    atom: usize,
    serial: u64,
    jitter: Seed,
    placement: Placement,
    buffer: Vec<u8>,
    cursor: usize,
    started: bool,
}

impl SyntheticCifSource {
    /// Streams `copies` copies of `tile`.
    #[must_use]
    pub fn new(tile: Tile, copies: u64) -> Self {
        let jitter = tile.jitter_seed(0);
        let placement = tile.placement(0, copies);
        Self {
            tile,
            copies,
            copy: 0,
            atom: 0,
            serial: 1,
            jitter,
            placement,
            buffer: Vec::new(),
            cursor: 0,
            started: false,
        }
    }

    /// Streams however many copies are needed to reach at least `atoms` atoms.
    #[must_use]
    pub fn with_atoms(tile: Tile, atoms: u64) -> Self {
        let copies = tile.copies_for(atoms);
        Self::new(tile, copies)
    }

    /// The number of coordinate rows this stream will emit.
    #[must_use]
    pub fn atom_count(&self) -> u64 {
        self.tile.atoms_in(self.copies)
    }

    /// How many copies of `tile` come to approximately `bytes` of mmCIF.
    ///
    /// The row width is measured by generating one copy rather than assumed,
    /// so the estimate follows the tile's own identifier lengths instead of a
    /// constant that would drift as the fixtures change.
    #[must_use]
    pub fn copies_for_bytes(tile: &Tile, bytes: u64) -> u64 {
        let mut probe = Self::new(tile.clone(), 1);
        let mut sink = vec![0_u8; 1 << 16];
        let mut measured = 0_u64;
        loop {
            match probe.read(&mut sink) {
                Ok(0) | Err(_) => break,
                Ok(count) => match u64::try_from(count) {
                    Ok(count) => measured = measured.saturating_add(count),
                    Err(_) => break,
                },
            }
        }
        if measured == 0 {
            return 0;
        }
        bytes.div_ceil(measured).max(1)
    }

    /// Refills the buffer, returning false once every row has been emitted.
    fn refill(&mut self) -> bool {
        self.buffer.clear();
        self.cursor = 0;

        if !self.started {
            self.started = true;
            self.buffer.extend_from_slice(HEADER.as_bytes());
        }

        let atoms = self.tile.atoms();
        if atoms.is_empty() {
            return !self.buffer.is_empty();
        }

        for _ in 0..ROWS_PER_REFILL {
            if self.copy >= self.copies {
                break;
            }
            let Some(atom) = atoms.get(self.atom) else {
                break;
            };
            let position = self.placement.apply(atom.position, &mut self.jitter);

            self.buffer.extend_from_slice(b"ATOM ");
            push_u64(&mut self.buffer, self.serial);
            self.buffer.push(b' ');
            self.buffer
                .extend_from_slice(atom.element.symbol().as_bytes());
            self.buffer.push(b' ');
            self.buffer
                .extend_from_slice(self.tile.name(atom.atom_name).as_bytes());
            self.buffer.extend_from_slice(b" . ");
            self.buffer
                .extend_from_slice(self.tile.name(atom.comp_id).as_bytes());
            self.buffer.push(b' ');
            self.buffer
                .extend_from_slice(self.tile.name(atom.asym_id).as_bytes());
            push_u64(&mut self.buffer, self.copy);
            self.buffer.extend_from_slice(b" 1 ");
            push_i32(&mut self.buffer, atom.seq_id);
            for axis in position {
                self.buffer.push(b' ');
                push_fixed(&mut self.buffer, f64::from(axis), 3);
            }
            self.buffer.extend_from_slice(b" 1.00 ");
            push_fixed(&mut self.buffer, b_factor_of(atom.seq_id), 2);
            self.buffer.extend_from_slice(b" 1\n");

            self.serial = self.serial.saturating_add(1);
            self.atom += 1;
            if self.atom >= atoms.len() {
                self.atom = 0;
                self.copy += 1;
                self.jitter = self.tile.jitter_seed(self.copy);
                self.placement = self.tile.placement(self.copy, self.copies);
            }
        }

        !self.buffer.is_empty()
    }
}

impl Read for SyntheticCifSource {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.cursor >= self.buffer.len() && !self.refill() {
            return Ok(0);
        }
        let Some(pending) = self.buffer.get(self.cursor..) else {
            return Ok(0);
        };
        let count = pending.len().min(out.len());
        let Some(source) = pending.get(..count) else {
            return Ok(0);
        };
        let Some(target) = out.get_mut(..count) else {
            return Ok(0);
        };
        target.copy_from_slice(source);
        self.cursor += count;
        Ok(count)
    }
}

/// A plausible temperature factor, varying along the chain rather than fixed.
fn b_factor_of(seq_id: i32) -> f64 {
    let along = f64::from(seq_id.rem_euclid(64));
    20.0 + along * 0.5
}

#[cfg(test)]
#[path = "emit_tests.rs"]
mod tests;
