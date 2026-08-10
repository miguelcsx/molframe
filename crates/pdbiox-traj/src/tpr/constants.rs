//! Versioned GROMACS topology function vocabulary.

pub(super) const DIMENSIONS: usize = 3;
pub(super) const FORMAT_SIZE_FIELD: i32 = 119;
pub(super) const REMOVE_THOLE_RFAC: i32 = 127;
pub(super) const FUNCTION_COUNT: i32 = 95;
pub(super) const SUPPORTED_VERSIONS: [i32; 16] = [
    58, 73, 83, 100, 103, 110, 112, 116, 119, 122, 127, 129, 133, 134, 137, 138,
];

pub(super) const F_BONDS: i32 = 0;
pub(super) const F_G96BONDS: i32 = 1;
pub(super) const F_MORSE: i32 = 2;
pub(super) const F_CUBICBONDS: i32 = 3;
pub(super) const F_CONNBONDS: i32 = 4;
pub(super) const F_HARMONIC: i32 = 5;
pub(super) const F_FENEBONDS: i32 = 6;
pub(super) const F_TABBONDS: i32 = 7;
pub(super) const F_TABBONDSNC: i32 = 8;
pub(super) const F_RESTRBONDS: i32 = 9;
pub(super) const F_ANGLES: i32 = 10;
pub(super) const F_G96ANGLES: i32 = 11;
pub(super) const F_RESTRANGLES: i32 = 12;
pub(super) const F_LINEAR_ANGLES: i32 = 13;
pub(super) const F_CROSS_BOND_BONDS: i32 = 14;
pub(super) const F_CROSS_BOND_ANGLES: i32 = 15;
pub(super) const F_UREY_BRADLEY: i32 = 16;
pub(super) const F_QUARTIC_ANGLES: i32 = 17;
pub(super) const F_TABANGLES: i32 = 18;
pub(super) const F_PDIHS: i32 = 19;
pub(super) const F_RBDIHS: i32 = 20;
pub(super) const F_RESTRDIHS: i32 = 21;
pub(super) const F_CBTDIHS: i32 = 22;
pub(super) const F_FOURDIHS: i32 = 23;
pub(super) const F_IDIHS: i32 = 24;
pub(super) const F_PIDIHS: i32 = 25;
pub(super) const F_TABDIHS: i32 = 26;
pub(super) const F_CMAP: i32 = 27;
pub(super) const F_GB12: i32 = 28;
pub(super) const F_GB13: i32 = 29;
pub(super) const F_GB14: i32 = 30;
pub(super) const F_LJ14: i32 = 33;
pub(super) const F_LJC14_Q: i32 = 35;
pub(super) const F_LJC_PAIRS_NB: i32 = 36;
pub(super) const F_LJ: i32 = 37;
pub(super) const F_BHAM: i32 = 38;
pub(super) const F_POLARIZATION: i32 = 48;
pub(super) const F_WATER_POL: i32 = 49;
pub(super) const F_THOLE_POL: i32 = 50;
pub(super) const F_ANHARM_POL: i32 = 51;
pub(super) const F_POSRES: i32 = 52;
pub(super) const F_FBPOSRES: i32 = 53;
pub(super) const F_DISRES: i32 = 54;
pub(super) const F_ORIRES: i32 = 56;
pub(super) const F_ANGRES: i32 = 58;
pub(super) const F_ANGRESZ: i32 = 59;
pub(super) const F_DIHRES: i32 = 60;
pub(super) const F_CONSTR: i32 = 62;
pub(super) const F_CONSTRNC: i32 = 63;
pub(super) const F_SETTLE: i32 = 64;
pub(super) const F_VSITE1: i32 = 65;
pub(super) const F_VSITE2: i32 = 66;
pub(super) const F_VSITE2FD: i32 = 67;
pub(super) const F_VSITE3: i32 = 68;
pub(super) const F_VSITE3FD: i32 = 69;
pub(super) const F_VSITE3FAD: i32 = 70;
pub(super) const F_VSITE3OUT: i32 = 71;
pub(super) const F_VSITE4FD: i32 = 72;
pub(super) const F_VSITE4FDN: i32 = 73;
pub(super) const F_VSITEN: i32 = 74;
pub(super) const F_ENNPOT: i32 = 78;

pub(super) const FUNCTION_ADDITIONS: [(i32, i32); 59] = [
    (20, F_CUBICBONDS),
    (20, F_CONNBONDS),
    (20, F_HARMONIC),
    (34, F_FENEBONDS),
    (43, F_TABBONDS),
    (43, F_TABBONDSNC),
    (70, F_RESTRBONDS),
    (98, F_RESTRANGLES),
    (76, F_LINEAR_ANGLES),
    (30, F_CROSS_BOND_BONDS),
    (30, F_CROSS_BOND_ANGLES),
    (30, F_UREY_BRADLEY),
    (34, F_QUARTIC_ANGLES),
    (43, F_TABANGLES),
    (98, F_RESTRDIHS),
    (98, F_CBTDIHS),
    (26, F_FOURDIHS),
    (26, F_PIDIHS),
    (43, F_TABDIHS),
    (65, F_CMAP),
    (60, F_GB12),
    (61, F_GB13),
    (61, F_GB14),
    (72, 31),
    (72, 32),
    (41, F_LJC14_Q),
    (41, F_LJC_PAIRS_NB),
    (32, 40),
    (32, 44),
    (32, 45),
    (93, 46),
    (46, 47),
    (30, F_POLARIZATION),
    (36, F_THOLE_POL),
    (90, F_FBPOSRES),
    (22, 55),
    (22, F_ORIRES),
    (22, 57),
    (26, F_DIHRES),
    (26, 61),
    (49, F_VSITE4FDN),
    (50, F_VSITEN),
    (46, 75),
    (20, 77),
    (46, 82),
    (69, 84),
    (66, 85),
    (54, 87),
    (76, F_ANHARM_POL),
    (79, 90),
    (79, 91),
    (79, 92),
    (79, 93),
    (79, 94),
    (117, 76),
    (121, F_VSITE1),
    (118, F_VSITE2FD),
    (137, F_ENNPOT),
    (79, 95),
];

pub(super) fn function_absent(version: i32, function: i32) -> bool {
    FUNCTION_ADDITIONS
        .iter()
        .any(|&(introduced, candidate)| version < introduced && function == candidate)
}

pub(super) fn remap_function(version: i32, mut function: i32) -> i32 {
    for &(introduced, candidate) in &FUNCTION_ADDITIONS {
        if version < introduced && function >= candidate {
            function += 1;
        }
    }
    function
}

pub(super) const fn interaction_arity(function: i32) -> Option<usize> {
    match function {
        0..=9 | 28..=30 | 33 | 35..=38 | 48 | 54 | 56 | 59 | 62 | 63 | 65 | 74 => Some(2),
        10..=18 | 64 | 66 | 67 => Some(3),
        19..=26 | 50 | 58 | 60 | 68..=71 => Some(4),
        27 | 49 | 72 | 73 => Some(5),
        _ => None,
    }
}

pub(super) const fn is_bond_function(function: i32) -> bool {
    matches!(function, 0..=9 | 62 | 63)
}
