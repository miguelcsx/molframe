//! Embedded structure and chemistry fixture definitions.

/// A structure size class, used to sweep a bench across realistic inputs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Sample {
    /// 1CRN crambin — 46 residues, ~327 atoms. Inner-loop signal.
    Tiny,
    /// 1UBQ ubiquitin — one chain plus ordered waters, ~660 atoms.
    Small,
    /// 4HHB haemoglobin — four chains and four HEM ligands, ~4.4k atoms.
    Medium,
    /// 1AON GroEL/GroES — 14 subunits, ~58k atoms. Scaling stress.
    Large,
    /// 2M7C ubiquitin — a 32-model NMR ensemble. Multi-model paths.
    Ensemble,
}

impl Sample {
    /// Every size class, smallest first.
    pub const ALL: [Sample; 5] = [
        Sample::Tiny,
        Sample::Small,
        Sample::Medium,
        Sample::Large,
        Sample::Ensemble,
    ];

    /// The single-model size classes, smallest first (excludes [`Sample::Ensemble`]).
    pub const MODELS: [Sample; 4] = [Sample::Tiny, Sample::Small, Sample::Medium, Sample::Large];

    /// A short stable label for Criterion parameter axes.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Sample::Tiny => "1crn",
            Sample::Small => "1ubq",
            Sample::Medium => "4hhb",
            Sample::Large => "1aon",
            Sample::Ensemble => "2m7c",
        }
    }

    /// The `BinaryCIF` bytes for this sample. Present for every size class.
    #[must_use]
    pub const fn bcif(self) -> &'static [u8] {
        match self {
            Sample::Tiny => include_bytes!("../data/1crn.bcif"),
            Sample::Small => include_bytes!("../data/1ubq.bcif"),
            Sample::Medium => include_bytes!("../data/4hhb.bcif"),
            Sample::Large => include_bytes!("../data/1aon.bcif"),
            Sample::Ensemble => include_bytes!("../data/2m7c.pdb"),
        }
    }

    /// The legacy PDB bytes, for the size classes that ship one.
    ///
    /// [`Sample::Large`] has none (committed as compact `BinaryCIF`).
    #[must_use]
    pub const fn pdb(self) -> Option<&'static [u8]> {
        match self {
            Sample::Tiny => Some(include_bytes!("../data/1crn.pdb")),
            Sample::Small => Some(include_bytes!("../data/1ubq.pdb")),
            Sample::Medium => Some(include_bytes!("../data/4hhb.pdb")),
            Sample::Ensemble => Some(include_bytes!("../data/2m7c.pdb")),
            Sample::Large => None,
        }
    }

    /// The uncompressed mmCIF bytes, for the size classes that ship one.
    ///
    /// [`Sample::Large`] ships mmCIF only gzip-compressed; see [`large_cif_gz`].
    /// The ensemble fixture ships only a PDB file.
    #[must_use]
    pub const fn cif(self) -> Option<&'static [u8]> {
        match self {
            Sample::Tiny => Some(include_bytes!("../data/1crn.cif")),
            Sample::Small => Some(include_bytes!("../data/1ubq.cif")),
            Sample::Medium => Some(include_bytes!("../data/4hhb.cif")),
            Sample::Ensemble | Sample::Large => None,
        }
    }
}

/// The gzip-compressed mmCIF for [`Sample::Large`].
#[must_use]
pub const fn large_cif_gz() -> &'static [u8] {
    include_bytes!("../data/1aon.cif.gz")
}

/// The Chemical Component Dictionary entry for haem (HEM).
#[must_use]
pub const fn ccd_hem() -> &'static [u8] {
    include_bytes!("../data/HEM.cif")
}

/// The Chemical Component Dictionary entry for ATP.
#[must_use]
pub const fn ccd_atp() -> &'static [u8] {
    include_bytes!("../data/ATP.cif")
}
