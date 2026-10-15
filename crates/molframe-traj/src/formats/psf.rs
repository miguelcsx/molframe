//! CHARMM/XPLOR/NAMD PSF topology reader and writer.

use std::fmt::Write;

/// One atom record from a PSF topology.
#[derive(Clone, Debug, PartialEq)]
pub struct PsfAtom {
    /// Segment identifier.
    pub segment: Box<str>,
    /// Residue identifier, retained as text for insertion-code variants.
    pub residue_id: Box<str>,
    /// Residue name.
    pub residue_name: Box<str>,
    /// Atom name.
    pub atom_name: Box<str>,
    /// Force-field atom type, textual or numeric.
    pub atom_type: Box<str>,
    /// Partial charge in proton-charge units.
    pub charge: f64,
    /// Mass in daltons.
    pub mass: f64,
}

/// PSF atoms and bonded interaction index tuples.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PsfTopology {
    /// Human-readable title records.
    pub titles: Vec<Box<str>>,
    /// Atom records in PSF order.
    pub atoms: Vec<PsfAtom>,
    /// Bond pairs.
    pub bonds: Vec<[u32; 2]>,
    /// Angle triples.
    pub angles: Vec<[u32; 3]>,
    /// Proper dihedral quadruples.
    pub dihedrals: Vec<[u32; 4]>,
    /// Improper dihedral quadruples.
    pub impropers: Vec<[u32; 4]>,
    /// Hydrogen-bond donor pairs.
    pub donors: Vec<[u32; 2]>,
    /// Hydrogen-bond acceptor pairs.
    pub acceptors: Vec<[u32; 2]>,
    /// CMAP cross-term octuples.
    pub cross_terms: Vec<[u32; 8]>,
}

/// Malformed or inconsistent PSF input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum PsfError {
    /// The PSF signature was absent.
    #[error("PSF signature is absent")]
    MissingSignature,
    /// A required section count was absent or invalid.
    #[error("invalid PSF section count")]
    InvalidSectionCount,
    /// An atom record was incomplete or non-numeric where required.
    #[error("invalid PSF atom record")]
    InvalidAtom,
    /// A connectivity section did not contain its declared number of indices.
    #[error("truncated PSF connectivity section")]
    TruncatedSection,
    /// A connectivity index is zero or beyond the declared atom count.
    #[error("PSF connectivity index is outside the atom table")]
    IndexOutOfRange,
}

/// Parses standard, `EXT`, XPLOR and NAMD text PSF files.
///
/// # Errors
///
/// Returns an error for invalid counts, atom records, truncated connectivity or
/// one-based indices outside `NATOM`.
pub fn parse_psf(text: &str) -> Result<PsfTopology, PsfError> {
    let lines: Vec<_> = text.lines().collect();
    if lines
        .first()
        .is_none_or(|line| line.split_whitespace().next() != Some("PSF"))
    {
        return Err(PsfError::MissingSignature);
    }
    let mut topology = PsfTopology::default();
    let mut cursor = 1;
    while cursor < lines.len() {
        let line = lines[cursor];
        cursor += 1;
        let Some(marker) = marker(line) else {
            continue;
        };
        let count = section_count(line)?;
        match marker {
            "NTITLE" => {
                for title in take_lines(&lines, &mut cursor, count)? {
                    let title = title.trim_start_matches('*').trim();
                    let title = match title.strip_prefix("REMARKS") {
                        Some(remark) => remark,
                        None => title,
                    };
                    topology.titles.push(title.trim().into());
                }
            }
            "NATOM" => {
                for atom in take_lines(&lines, &mut cursor, count)? {
                    topology.atoms.push(parse_atom(atom)?);
                }
            }
            "NBOND" => topology.bonds = tuples(&lines, &mut cursor, count, topology.atoms.len())?,
            "NTHETA" => {
                topology.angles = tuples(&lines, &mut cursor, count, topology.atoms.len())?;
            }
            "NPHI" => {
                topology.dihedrals = tuples(&lines, &mut cursor, count, topology.atoms.len())?;
            }
            "NIMPHI" => {
                topology.impropers = tuples(&lines, &mut cursor, count, topology.atoms.len())?;
            }
            "NDON" => topology.donors = tuples(&lines, &mut cursor, count, topology.atoms.len())?,
            "NACC" => {
                topology.acceptors = tuples(&lines, &mut cursor, count, topology.atoms.len())?;
            }
            "NCRTERM" => {
                topology.cross_terms = tuples(&lines, &mut cursor, count, topology.atoms.len())?;
            }
            _ => {}
        }
    }
    if topology.atoms.is_empty() {
        return Err(PsfError::InvalidSectionCount);
    }
    Ok(topology)
}

fn marker(line: &str) -> Option<&str> {
    line.split_once('!')?
        .1
        .split_whitespace()
        .next()
        .map(|value| value.trim_end_matches(':'))
}

fn section_count(line: &str) -> Result<usize, PsfError> {
    line.split_whitespace()
        .next()
        .and_then(|value| value.parse().ok())
        .ok_or(PsfError::InvalidSectionCount)
}

fn take_lines<'a>(
    lines: &'a [&str],
    cursor: &mut usize,
    count: usize,
) -> Result<&'a [&'a str], PsfError> {
    let end = cursor
        .checked_add(count)
        .filter(|end| *end <= lines.len())
        .ok_or(PsfError::TruncatedSection)?;
    let records = &lines[*cursor..end];
    *cursor = end;
    Ok(records)
}

fn parse_atom(line: &str) -> Result<PsfAtom, PsfError> {
    let fields: Vec<_> = line.split_whitespace().collect();
    if fields.len() < 8 || fields[0].parse::<u32>().ok().is_none_or(|value| value == 0) {
        return Err(PsfError::InvalidAtom);
    }
    Ok(PsfAtom {
        segment: fields[1].into(),
        residue_id: fields[2].into(),
        residue_name: fields[3].into(),
        atom_name: fields[4].into(),
        atom_type: fields[5].into(),
        charge: fields[6].parse().map_err(|_| PsfError::InvalidAtom)?,
        mass: fields[7].parse().map_err(|_| PsfError::InvalidAtom)?,
    })
}

fn tuples<const N: usize>(
    lines: &[&str],
    cursor: &mut usize,
    count: usize,
    atom_count: usize,
) -> Result<Vec<[u32; N]>, PsfError> {
    let needed = count.checked_mul(N).ok_or(PsfError::InvalidSectionCount)?;
    let mut values = Vec::with_capacity(needed);
    while values.len() < needed && *cursor < lines.len() && marker(lines[*cursor]).is_none() {
        for value in lines[*cursor].split_whitespace() {
            let index: u32 = value.parse().map_err(|_| PsfError::TruncatedSection)?;
            if index == 0
                || usize::try_from(index).map_or(true, |atom_index| atom_index > atom_count)
            {
                return Err(PsfError::IndexOutOfRange);
            }
            values.push(index - 1);
        }
        *cursor += 1;
    }
    if values.len() != needed {
        return Err(PsfError::TruncatedSection);
    }
    let (rows, _) = values.as_chunks::<N>();
    Ok(rows
        .iter()
        .map(|chunk| std::array::from_fn(|index| chunk[index]))
        .collect())
}

/// Writes an XPLOR-compatible PSF with deterministic section order.
#[must_use]
pub fn write_psf(topology: &PsfTopology) -> String {
    let mut output = String::from("PSF XPLOR\n\n");
    let _ = writeln!(output, "{:>8} !NTITLE", topology.titles.len());
    for title in &topology.titles {
        let _ = writeln!(output, " REMARKS {title}");
    }
    let _ = writeln!(output, "\n{:>8} !NATOM", topology.atoms.len());
    for (index, atom) in topology.atoms.iter().enumerate() {
        let _ = writeln!(
            output,
            "{:>8} {:<8} {:<8} {:<8} {:<8} {:<8} {:>14.6} {:>13.4}           0",
            index + 1,
            atom.segment,
            atom.residue_id,
            atom.residue_name,
            atom.atom_name,
            atom.atom_type,
            atom.charge,
            atom.mass,
        );
    }
    write_tuples(&mut output, "NBOND: bonds", &topology.bonds);
    write_tuples(&mut output, "NTHETA: angles", &topology.angles);
    write_tuples(&mut output, "NPHI: dihedrals", &topology.dihedrals);
    write_tuples(&mut output, "NIMPHI: impropers", &topology.impropers);
    write_tuples(&mut output, "NDON: donors", &topology.donors);
    write_tuples(&mut output, "NACC: acceptors", &topology.acceptors);
    write_tuples(&mut output, "NCRTERM: cross-terms", &topology.cross_terms);
    output
}

fn write_tuples<const N: usize>(output: &mut String, marker: &str, tuples: &[[u32; N]]) {
    let _ = writeln!(output, "\n{:>8} !{marker}", tuples.len());
    for tuple in tuples {
        for index in tuple {
            let _ = write!(output, "{:>8}", index + 1);
        }
        output.push('\n');
    }
}

#[cfg(test)]
#[path = "psf_tests.rs"]
mod tests;
