//! Version-aware force-field parameter striding.

use super::constants::{
    F_ANGLES, F_ANGRES, F_ANGRESZ, F_ANHARM_POL, F_BHAM, F_BONDS, F_CBTDIHS, F_CMAP, F_CONNBONDS,
    F_CONSTR, F_CONSTRNC, F_CROSS_BOND_ANGLES, F_CROSS_BOND_BONDS, F_CUBICBONDS, F_DIHRES,
    F_DISRES, F_ENNPOT, F_FBPOSRES, F_FENEBONDS, F_FOURDIHS, F_G96ANGLES, F_G96BONDS, F_GB12,
    F_GB13, F_GB14, F_HARMONIC, F_IDIHS, F_LINEAR_ANGLES, F_LJ, F_LJ14, F_LJC_PAIRS_NB, F_LJC14_Q,
    F_MORSE, F_ORIRES, F_PDIHS, F_PIDIHS, F_POLARIZATION, F_POSRES, F_QUARTIC_ANGLES, F_RBDIHS,
    F_RESTRANGLES, F_RESTRBONDS, F_RESTRDIHS, F_SETTLE, F_TABANGLES, F_TABBONDS, F_TABBONDSNC,
    F_TABDIHS, F_THOLE_POL, F_UREY_BRADLEY, F_VSITE1, F_VSITE2, F_VSITE2FD, F_VSITE3, F_VSITE3FAD,
    F_VSITE3FD, F_VSITE3OUT, F_VSITE4FD, F_VSITE4FDN, F_VSITEN, F_WATER_POL, REMOVE_THOLE_RFAC,
};
use super::error::TprError;
use super::xdr::Decoder;

pub(super) fn skip_parameter_table(
    decoder: &mut Decoder<'_>,
    functions: &[i32],
    version: i32,
) -> Result<(), TprError> {
    for &function in functions {
        skip_parameter(decoder, function, version)?;
    }
    Ok(())
}

fn skip_parameter(decoder: &mut Decoder<'_>, function: i32, version: i32) -> Result<(), TprError> {
    match function {
        F_ANGLES | F_G96ANGLES | F_BONDS | F_G96BONDS | F_HARMONIC | F_IDIHS => {
            decoder.skip_reals(4)
        }
        F_RESTRANGLES | F_RESTRDIHS => skip_versioned_reals(decoder, version, 2, 134, 2),
        F_LINEAR_ANGLES | F_CROSS_BOND_ANGLES | F_LJ14 | F_LJC_PAIRS_NB => decoder.skip_reals(4),
        F_FENEBONDS | F_LJ | F_CONSTR | F_CONSTRNC | F_SETTLE | F_VSITE3 | F_VSITE3FD
        | F_VSITE3FAD => decoder.skip_reals(2),
        F_RESTRBONDS => decoder.skip_reals(8),
        F_TABBONDS | F_TABBONDSNC | F_TABANGLES | F_TABDIHS => {
            decoder.skip_reals(1)?;
            decoder.i32()?;
            decoder.skip_reals(1)
        }
        F_CROSS_BOND_BONDS | F_BHAM | F_CUBICBONDS | F_ANHARM_POL | F_VSITE3OUT | F_VSITE4FD
        | F_VSITE4FDN => decoder.skip_reals(3),
        F_UREY_BRADLEY => skip_versioned_reals(decoder, version, 4, 79, 4),
        F_QUARTIC_ANGLES | F_WATER_POL => decoder.skip_reals(6),
        F_MORSE => skip_versioned_reals(decoder, version, 3, 79, 3),
        F_CONNBONDS | F_VSITE1 | F_ENNPOT => Ok(()),
        F_POLARIZATION | F_VSITE2 | F_VSITE2FD => decoder.skip_reals(1),
        F_THOLE_POL => decoder.skip_reals(if version < REMOVE_THOLE_RFAC { 4 } else { 3 }),
        F_LJC14_Q => decoder.skip_reals(5),
        F_PIDIHS | F_ANGRES | F_ANGRESZ | F_PDIHS => {
            decoder.skip_reals(4)?;
            decoder.i32().map(|_| ())
        }
        F_DISRES => {
            decoder.i32()?;
            decoder.i32()?;
            decoder.skip_reals(4)
        }
        F_ORIRES => {
            skip_ints(decoder, 3)?;
            decoder.skip_reals(3)
        }
        F_DIHRES => skip_dihedral_restraint(decoder, version),
        F_POSRES | F_RBDIHS | F_FOURDIHS => decoder.skip_reals(12),
        F_FBPOSRES => {
            decoder.i32()?;
            decoder.skip_reals(5)
        }
        F_CBTDIHS => skip_versioned_reals(decoder, version, 6, 134, 6),
        F_VSITEN => {
            decoder.i32()?;
            decoder.skip_reals(1)
        }
        F_GB12 | F_GB13 | F_GB14 => {
            if version < 68 {
                decoder.skip_reals(4)?;
            }
            decoder.skip_reals(5)
        }
        F_CMAP => {
            decoder.i32()?;
            decoder.i32().map(|_| ())
        }
        _ => Err(TprError::UnsupportedFunction(function)),
    }
}

fn skip_versioned_reals(
    decoder: &mut Decoder<'_>,
    version: i32,
    base: usize,
    threshold: i32,
    added: usize,
) -> Result<(), TprError> {
    decoder.skip_reals(base)?;
    if version >= threshold {
        decoder.skip_reals(added)?;
    }
    Ok(())
}

fn skip_dihedral_restraint(decoder: &mut Decoder<'_>, version: i32) -> Result<(), TprError> {
    if version < 72 {
        skip_ints(decoder, 2)?;
    }
    decoder.skip_reals(3)?;
    if version >= 72 {
        decoder.skip_reals(3)?;
    }
    Ok(())
}

fn skip_ints(decoder: &mut Decoder<'_>, count: usize) -> Result<(), TprError> {
    for _ in 0..count {
        decoder.i32()?;
    }
    Ok(())
}
