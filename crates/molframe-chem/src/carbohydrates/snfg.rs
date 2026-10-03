//! Curated Mol* SNFG metadata and common CCD aliases, adapted under MIT.
//!
//! Copyright (c) 2018–2026 Mol* contributors. The retained permission notice is
//! in data/LICENSE.molstar; data/README.md records the pinned source and scope.

use super::types::{SnfgShape, SnfgSymbol};

const BLUE: u32 = 0x00_90_bc;
const GREEN: u32 = 0x00_a6_51;
const YELLOW: u32 = 0xff_d4_00;
const ORANGE: u32 = 0xf4_79_20;
const PINK: u32 = 0xf6_9e_a1;
const PURPLE: u32 = 0xa5_43_99;
const LIGHT_BLUE: u32 = 0x8f_cc_e9;
const BROWN: u32 = 0xa1_7a_4d;
const RED: u32 = 0xed_1c_24;
const SECONDARY: u32 = 0xf1_ec_e1;

pub(super) const UNKNOWN: SnfgSymbol = SnfgSymbol {
    abbreviation: "Unk",
    name: "Unknown",
    color: SECONDARY,
    shape: SnfgShape::FlatHexagon,
    secondary_color: None,
};

/// Looks up a common CCD identifier or case-sensitive SNFG abbreviation.
///
/// The curated CCD aliases exclude force-field residue names, whose namespaces
/// collide with CCD identifiers. Unknown identifiers return `None`; topology
/// alone cannot distinguish glucose from a stereoisomer.
///
/// ```
/// use molframe_chem::{snfg_symbol, SnfgShape};
/// assert_eq!(snfg_symbol("NAG").map(|symbol| symbol.abbreviation), Some("GlcNAc"));
/// assert_eq!(snfg_symbol("NAG").map(|symbol| symbol.shape), Some(SnfgShape::FilledCube));
/// ```
#[must_use]
pub fn snfg_symbol(component: &str) -> Option<SnfgSymbol> {
    let abbreviation = match component {
        "GLC" | "BGC" | "Z8T" | "TRE" | "MLR" => "Glc",
        "MAN" | "BMA" => "Man",
        "GLA" | "GAL" | "GZL" | "GXL" | "GIV" => "Gal",
        "4GL" | "GL0" | "GUP" | "Z8H" => "Gul",
        "Z6H" | "3MK" | "SHD" => "Alt",
        "AFD" | "ALL" | "WOO" | "Z2D" => "All",
        "ZEE" | "A5C" | "SDY" => "Tal",
        "ZCD" | "Z0F" | "4N2" => "Ido",
        "NDG" | "NAG" | "NGZ" => "GlcNAc",
        "BM3" | "BM7" => "ManNAc",
        "A2G" | "NGA" | "YYQ" => "GalNAc",
        "LXB" => "GulNAc",
        "NAA" => "AllNAc",
        "LXZ" | "HSQ" => "IdoNAc",
        "PA1" | "GCS" => "GlcN",
        "95Z" => "ManN",
        "X6X" | "1GN" => "GalN",
        "GCU" | "BDP" => "GlcA",
        "MAV" | "BEM" => "ManA",
        "ADA" | "GTR" | "GTK" => "GalA",
        "LGU" => "GulA",
        "X1X" | "X0X" => "TalA",
        "IDR" => "IdoA",
        "G6D" | "YYK" => "Qui",
        "RAM" | "RM4" | "XXR" => "Rha",
        "66O" => "6dGul",
        "FUC" | "FUL" | "FCA" | "FCB" | "GYE" => "Fuc",
        "Z9W" => "QuiNAc",
        "49T" => "FucNAc",
        "DDA" | "RAE" | "Z5J" => "Oli",
        "TYV" => "Tyv",
        "ABE" => "Abe",
        "PZU" => "Par",
        "Z3U" => "Dig",
        "64K" | "ARA" | "ARB" | "AHR" | "FUB" | "BXY" | "BXX" | "SEJ" => "Ara",
        "LDY" | "Z4W" => "Lyx",
        "XYS" | "XYP" | "XYZ" | "HSY" | "LXC" => "Xyl",
        "YYM" | "RIP" | "RIB" | "BDR" | "0MK" | "Z6J" | "32O" => "Rib",
        "KDM" | "KDN" => "Kdn",
        "SIA" | "SLB" => "Neu5Ac",
        "NGC" | "NGE" => "Neu5Gc",
        "GMH" => "LDmanHep",
        "KDO" => "Kdo",
        "289" => "DDmanHep",
        "MUB" | "AMU" => "MurNAc",
        "1S4" | "MUR" => "Mur",
        "XXM" => "Api",
        "BDF" | "Z9N" | "FRU" | "LFR" => "Fru",
        "T6T" => "Tag",
        "SOE" | "UEA" => "Sor",
        "PSV" | "SF6" | "SF9" | "TTV" => "Psi",
        other => other,
    };
    symbol_for_abbreviation(abbreviation)
}

use SnfgShape::{
    CrossedCube, DividedCone, DividedDiamond, FilledCone, FilledCube, FilledDiamond, FilledSphere,
    FilledStar, FlatBox, FlatDiamond, FlatHexagon, Pentagon,
};

const SYMBOLS: &[(&str, &str, u32, SnfgShape)] = &[
    ("Glc", "Glucose", BLUE, FilledSphere),
    ("Man", "Mannose", GREEN, FilledSphere),
    ("Gal", "Galactose", YELLOW, FilledSphere),
    ("Gul", "Gulose", ORANGE, FilledSphere),
    ("Alt", "Altrose", PINK, FilledSphere),
    ("All", "Allose", PURPLE, FilledSphere),
    ("Tal", "Talose", LIGHT_BLUE, FilledSphere),
    ("Ido", "Idose", BROWN, FilledSphere),
    ("GlcNAc", "N-Acetyl Glucosamine", BLUE, FilledCube),
    ("ManNAc", "N-Acetyl Mannosamine", GREEN, FilledCube),
    ("GalNAc", "N-Acetyl Galactosamine", YELLOW, FilledCube),
    ("GulNAc", "N-Acetyl Gulosamine", ORANGE, FilledCube),
    ("AltNAc", "N-Acetyl Altrosamine", PINK, FilledCube),
    ("AllNAc", "N-Acetyl Allosamine", PURPLE, FilledCube),
    ("TalNAc", "N-Acetyl Talosamine", LIGHT_BLUE, FilledCube),
    ("IdoNAc", "N-Acetyl Idosamine", BROWN, FilledCube),
    ("GlcN", "Glucosamine", BLUE, CrossedCube),
    ("ManN", "Mannosamine", GREEN, CrossedCube),
    ("GalN", "Galactosamine", YELLOW, CrossedCube),
    ("GulN", "Gulosamine", ORANGE, CrossedCube),
    ("AltN", "Altrosamine", PINK, CrossedCube),
    ("AllN", "Allosamine", PURPLE, CrossedCube),
    ("TalN", "Talosamine", LIGHT_BLUE, CrossedCube),
    ("IdoN", "Idosamine", BROWN, CrossedCube),
    ("GlcNS", "N-sulfo Glucosamine", BLUE, CrossedCube),
    ("GlcA", "Glucuronic Acid", BLUE, DividedDiamond),
    ("ManA", "Mannuronic Acid", GREEN, DividedDiamond),
    ("GalA", "Galacturonic Acid", YELLOW, DividedDiamond),
    ("GulA", "Guluronic Acid", ORANGE, DividedDiamond),
    ("AltA", "Altruronic Acid", PINK, DividedDiamond),
    ("AllA", "Alluronic Acid", PURPLE, DividedDiamond),
    ("TalA", "Taluronic Acid", LIGHT_BLUE, DividedDiamond),
    ("IdoA", "Iduronic Acid", BROWN, DividedDiamond),
    ("Qui", "Quinovose", BLUE, FilledCone),
    ("Rha", "Rhamnose", GREEN, FilledCone),
    ("6dGul", "6-Deoxy Gulose", ORANGE, FilledCone),
    ("6dAlt", "6-Deoxy Altrose", PINK, FilledCone),
    ("6dTal", "6-Deoxy Talose", LIGHT_BLUE, FilledCone),
    ("Fuc", "Fucose", RED, FilledCone),
    ("QuiNAc", "N-Acetyl Quinovosamine", BLUE, DividedCone),
    ("RhaNAc", "N-Acetyl Rhamnosamine", GREEN, DividedCone),
    (
        "6dAltNAc",
        "N-Acetyl 6-Deoxy Altrosamine",
        PINK,
        DividedCone,
    ),
    (
        "6dTalNAc",
        "N-Acetyl 6-Deoxy Talosamine",
        LIGHT_BLUE,
        DividedCone,
    ),
    ("FucNAc", "N-Acetyl Fucosamine", RED, DividedCone),
    ("AAT", "AAT", RED, DividedCone),
    ("Oli", "Olivose", BLUE, FlatBox),
    ("Tyv", "Tyvelose", GREEN, FlatBox),
    ("Abe", "Abequose", ORANGE, FlatBox),
    ("Par", "Paratose", PINK, FlatBox),
    ("Dig", "Digitoxose", PURPLE, FlatBox),
    ("Col", "Colitose", LIGHT_BLUE, FlatBox),
    ("Ara", "Arabinose", GREEN, FilledStar),
    ("Lyx", "Lyxose", YELLOW, FilledStar),
    ("Xyl", "Xylose", ORANGE, FilledStar),
    ("Rib", "Ribose", PINK, FilledStar),
    ("Kdn", "Keto-Deoxy Nonulonic Acid", GREEN, FilledDiamond),
    ("Neu5Ac", "N-Acetyl Neuraminic Acid", PURPLE, FilledDiamond),
    (
        "Neu5Gc",
        "N-Glycolyl Neuraminic Acid",
        LIGHT_BLUE,
        FilledDiamond,
    ),
    ("Neu", "Neuraminic Acid", BROWN, FilledDiamond),
    ("Sia", "Sialic acid", RED, FilledDiamond),
    ("Pse", "Pseudaminic Acid", GREEN, FlatDiamond),
    ("Leg", "Legionaminic Acid", YELLOW, FlatDiamond),
    ("Aci", "Acinetaminic Acid", PINK, FlatDiamond),
    ("4eLeg", "4-Epilegionaminic Acid", LIGHT_BLUE, FlatDiamond),
    ("Bac", "Bacillosamine", BLUE, FlatHexagon),
    ("LDmanHep", "L-Glycero-D-Manno Heptose", GREEN, FlatHexagon),
    ("Kdo", "Keto-Deoxy Octulonic Acid", YELLOW, FlatHexagon),
    ("Dha", "3-Deoxy Lyxo-Heptulosaric Acid", ORANGE, FlatHexagon),
    ("DDmanHep", "D-Glycero-D-Manno-Heptose", PINK, FlatHexagon),
    ("MurNAc", "N-Acetyl Muramic Acid", PURPLE, FlatHexagon),
    ("MurNGc", "N-Glycolyl Muramic Acid", LIGHT_BLUE, FlatHexagon),
    ("Mur", "Muramic Acid", BROWN, FlatHexagon),
    (
        "K3O",
        "D-glycero-a-D-talo-oct-2-Ulosonic Acid",
        RED,
        FlatHexagon,
    ),
    (
        "dUA",
        "4-deoxy-4,5-didehydro iduronic acid",
        SECONDARY,
        FlatHexagon,
    ),
    ("Api", "Apicose", GREEN, Pentagon),
    ("Fru", "Fructose", GREEN, Pentagon),
    ("Tag", "Tagatose", YELLOW, Pentagon),
    ("Sor", "Sorbose", ORANGE, Pentagon),
    ("Psi", "Psicose", PINK, Pentagon),
];

fn symbol_for_abbreviation(value: &str) -> Option<SnfgSymbol> {
    let &(abbreviation, name, color, shape) = SYMBOLS
        .iter()
        .find(|(abbreviation, _, _, _)| *abbreviation == value)?;
    let secondary_color = match shape {
        CrossedCube | DividedDiamond | DividedCone => Some(SECONDARY),
        _ => None,
    };
    Some(SnfgSymbol {
        abbreviation,
        name,
        color,
        shape,
        secondary_color,
    })
}
