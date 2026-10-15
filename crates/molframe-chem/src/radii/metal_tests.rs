use super::*;
use crate::test_data::shipped;

const HEADER: &str = "name: lab\nversion: 1\nsource: survey\nlicence: CC0-1.0\n";

fn sites(body: &str) -> Provided<MetalSites> {
    match MetalSites::decode(
        format!("{HEADER}{body}").as_bytes(),
        TableIdentity::supplied("x", "x", "x", "x"),
    ) {
        Ok(provided) => provided,
        Err(error) => panic!("sites refused: {error}"),
    }
}

fn ionic(body: &str) -> Provided<IonicRadii> {
    match IonicRadii::decode(
        format!("{HEADER}{body}").as_bytes(),
        TableIdentity::supplied("x", "x", "x", "x"),
    ) {
        Ok(provided) => provided,
        Err(error) => panic!("ionic refused: {error}"),
    }
}

#[test]
fn a_site_gives_its_vdw_radius_and_the_ionic_radius_of_exactly_that_state() {
    let metals = sites("Fe 2 VI HS 2.44\n");
    let catalogue =
        ionic("Fe 2 VI LS 0.61 0.75 -\nFe 2 VI HS 0.78 0.92 yes\nFe 3 VI HS 0.645 0.785 yes\n");
    let table = &metals.table;
    assert_eq!(
        table.radius(Element::IRON, MetalRadiusKind::Vdw, None),
        Some(2.44)
    );
    assert_eq!(
        table.radius(
            Element::IRON,
            MetalRadiusKind::Ionic,
            Some(&catalogue.table)
        ),
        Some(0.78)
    );
    assert_eq!(
        table.radius(Element::IRON, MetalRadiusKind::Ionic, None),
        None
    );
    assert_eq!(
        table.radius(Element::CARBON, MetalRadiusKind::Vdw, None),
        None
    );
    assert!(table.is_metal(Element::IRON) && !table.is_metal(Element::CARBON));
    assert_eq!(table.elements().collect::<Vec<_>>(), vec![Element::IRON]);
}

#[test]
fn defects_are_refused() {
    for body in [
        "",
        "Fe 2 VI HS\n",
        "Fe 2 VI HS 0\n",
        "Fe 2 VI x 2.4\n",
        "Fe 2 VI - 2.4\nFe 3 VI - 2.4\n",
    ] {
        let refused = MetalSites::decode(
            format!("{HEADER}{body}").as_bytes(),
            TableIdentity::supplied("x", "x", "x", "x"),
        );
        assert!(refused.is_err(), "{body:?}");
    }
}

#[test]
fn the_shipped_typical_sites_exist_in_the_shipped_shannon_catalogue() {
    let metals = shipped::<MetalSites>("macromolecular-metals.metals");
    let catalogue = shipped::<IonicRadii>("shannon-ionic.ionic");
    for element in metals.table.elements() {
        let ionic = metals
            .table
            .radius(element, MetalRadiusKind::Ionic, Some(&catalogue.table));
        let vdw = metals.table.radius(element, MetalRadiusKind::Vdw, None);
        let (Some(ionic), Some(vdw)) = (ionic, vdw) else {
            panic!("{} lacks a radius", element.symbol())
        };
        assert!(
            ionic < vdw,
            "{}: a cation is smaller than its vdW sphere",
            element.symbol()
        );
    }
    assert_eq!(
        metals.table.radius(
            Element::IRON,
            MetalRadiusKind::Ionic,
            Some(&catalogue.table)
        ),
        Some(0.78)
    );
}
