use super::{METAL_SITES, MetalRadiusKind, metal_radius, metal_site};
use molframe_core::Element;

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).unwrap_or_else(|| panic!("unknown symbol {symbol}"))
}

fn ionic(symbol: &str) -> f32 {
    metal_radius(element(symbol), MetalRadiusKind::Ionic)
        .unwrap_or_else(|| panic!("{symbol} has no ionic radius"))
}

#[test]
fn ionic_radii_match_shannon_effective_values() {
    // Shannon (1976) effective radii at the typical site, in angstrom.
    for (symbol, expected) in [
        ("Na", 1.02),
        ("K", 1.38),
        ("Mg", 0.72),
        ("Ca", 1.00),
        ("Mn", 0.83),
        ("Fe", 0.78),
        ("Co", 0.745),
        ("Ni", 0.69),
        ("Cu", 0.65),
        ("Zn", 0.60),
    ] {
        assert!((ionic(symbol) - expected).abs() < 1e-4, "{symbol}");
    }
}

#[test]
fn van_der_waals_radii_are_alvarez() {
    for (symbol, expected) in [("Mg", 2.51), ("Ca", 2.62), ("Zn", 2.39), ("Fe", 2.44)] {
        let radius = metal_radius(element(symbol), MetalRadiusKind::Vdw);
        assert_eq!(radius, Some(expected), "{symbol}");
    }
}

#[test]
fn every_site_resolves_all_three_radii_in_order() {
    for site in &METAL_SITES {
        let e = element(site.symbol);
        let vdw = metal_radius(e, MetalRadiusKind::Vdw);
        let covalent = metal_radius(e, MetalRadiusKind::Covalent);
        let ionic = metal_radius(e, MetalRadiusKind::Ionic);
        let (Some(vdw), Some(covalent), Some(ionic)) = (vdw, covalent, ionic) else {
            panic!("{} lacks a radius", site.symbol)
        };
        assert!(
            ionic < vdw,
            "{}: cation is smaller than its vdW sphere",
            site.symbol
        );
        assert!(covalent < vdw, "{}", site.symbol);
    }
}

#[test]
fn bonding_reach_admits_the_ionic_contact_distance_to_oxide() {
    // The coordination window must contain the crystal contact distance, the
    // sum of the cation radius and Shannon's oxide radius (1.40 A), within
    // sampling tolerance, and stay below the van der Waals contact.
    for site in &METAL_SITES {
        let e = element(site.symbol);
        let contact = ionic(site.symbol) + 1.40;
        let window = (site.reach + 1.52) / 1.95;
        assert!(window >= contact - 0.05, "{} window {window}", site.symbol);
        let vdw = metal_radius(e, MetalRadiusKind::Vdw).unwrap_or(0.0);
        assert!(window < vdw + 1.52, "{}", site.symbol);
    }
}

#[test]
fn non_metals_have_no_site() {
    assert!(metal_site("C").is_none());
    assert_eq!(metal_radius(element("C"), MetalRadiusKind::Vdw), None);
}
