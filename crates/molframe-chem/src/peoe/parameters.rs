//! Versioned orbital-electronegativity parameter data.

use super::{PeoeAtomType, PeoeParameterProfile};

#[derive(Clone, Copy)]
pub(super) struct Coefficients {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl Coefficients {
    pub const fn ionisation_electronegativity(self, atom_type: PeoeAtomType) -> f64 {
        match atom_type {
            PeoeAtomType::H => 20.02,
            _ => self.a + self.b + self.c,
        }
    }
}

pub(super) const fn coefficients(
    profile: PeoeParameterProfile,
    atom_type: PeoeAtomType,
) -> Coefficients {
    match profile {
        PeoeParameterProfile::GasteigerMarsili => gasteiger_marsili(atom_type),
    }
}

const fn gasteiger_marsili(atom_type: PeoeAtomType) -> Coefficients {
    let values = match atom_type {
        PeoeAtomType::H => (7.17, 6.24, -0.56),
        PeoeAtomType::CSp3 => (7.98, 9.18, 1.88),
        PeoeAtomType::CSp2 => (8.79, 9.32, 1.51),
        PeoeAtomType::CSp => (10.39, 9.45, 0.73),
        PeoeAtomType::NSp3 => (11.54, 10.82, 1.36),
        PeoeAtomType::NSp2 => (12.87, 11.15, 0.85),
        PeoeAtomType::NSp => (15.68, 11.70, -0.27),
        PeoeAtomType::OSp3 => (14.18, 12.92, 1.39),
        PeoeAtomType::OSp2 => (17.07, 13.79, 0.47),
        PeoeAtomType::FSp3 => (14.66, 13.85, 2.31),
        PeoeAtomType::ClSp3 => (11.00, 9.69, 1.35),
        PeoeAtomType::BrSp3 => (10.08, 8.47, 1.16),
        PeoeAtomType::ISp3 => (9.90, 7.96, 0.96),
        PeoeAtomType::SSp3 | PeoeAtomType::SO => (10.14, 9.13, 1.38),
        PeoeAtomType::SO2 => (12.00, 10.81, 1.20),
        PeoeAtomType::SSp2 => (10.88, 9.49, 1.33),
        PeoeAtomType::PSp3 => (8.90, 8.24, 0.96),
        PeoeAtomType::PSp2 => (9.665, 8.530, 0.735),
        PeoeAtomType::SiSp3 => (7.300, 6.567, 0.657),
        PeoeAtomType::SiSp2 => (7.905, 6.748, 0.443),
        PeoeAtomType::SiSp => (9.065, 7.027, -0.002),
        PeoeAtomType::BSp3 => (5.980, 6.820, 1.605),
        PeoeAtomType::BSp2 => (6.420, 6.807, 1.322),
        PeoeAtomType::BeSp3 => (3.845, 6.755, 3.165),
        PeoeAtomType::BeSp2 => (4.005, 6.725, 3.035),
        PeoeAtomType::MgSp3 => (3.300, 5.587, 2.447),
        PeoeAtomType::MgSp2 => (3.565, 5.572, 2.197),
        PeoeAtomType::MgSp => (4.040, 5.472, 1.823),
        PeoeAtomType::AlSp3 => (5.375, 4.953, 0.867),
        PeoeAtomType::AlSp2 => (5.795, 5.020, 0.695),
    };
    Coefficients {
        a: values.0,
        b: values.1,
        c: values.2,
    }
}
