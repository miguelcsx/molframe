use super::*;

fn ring(elements: &[Element], orders: &[BondOrder], charges: &[Option<i8>]) -> Aromaticity {
    let bonds: Vec<_> = orders
        .iter()
        .enumerate()
        .map(|(i, &order)| (i, (i + 1) % elements.len(), order))
        .collect();
    perceive(elements, charges, &bonds)
}

#[test]
fn benzene_is_aromatic_in_both_encodings() {
    for orders in [
        vec![BondOrder::Aromatic; 6],
        vec![
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    ] {
        let result = ring(&[Element::CARBON; 6], &orders, &[Some(0); 6]);
        assert_eq!(result.atoms, vec![true; 6]);
        assert_eq!(result.bonds, vec![true; 6]);
    }
}

#[test]
fn pyrrole_and_furan_donate_a_lone_pair() {
    for hetero in [Element::NITROGEN, Element::OXYGEN, Element::SULFUR] {
        let result = ring(
            &[
                hetero,
                Element::CARBON,
                Element::CARBON,
                Element::CARBON,
                Element::CARBON,
            ],
            &[
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
            &[Some(0); 5],
        );
        assert!(result.atoms.iter().all(|value| *value));
    }
}

#[test]
fn saturated_interrupted_and_four_n_circuits_are_not_promoted() {
    for orders in [
        vec![BondOrder::Single; 6],
        vec![
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
        vec![
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    ] {
        let result = ring(
            &vec![Element::CARBON; orders.len()],
            &orders,
            &vec![Some(0); orders.len()],
        );
        assert!(result.atoms.iter().all(|value| !value));
    }
}

#[test]
fn exocyclic_double_bonds_and_unknown_charge_do_not_invent_aromaticity() {
    let mut bonds: Vec<_> = (0..6)
        .map(|i| {
            (
                i,
                (i + 1) % 6,
                if i % 2 == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                },
            )
        })
        .collect();
    bonds.push((0, 6, BondOrder::Double));
    assert!(!perceive(&[Element::CARBON; 7], &[Some(0); 7], &bonds).atoms[0]);
    bonds.pop();
    assert!(
        perceive(&[Element::CARBON; 6], &[None; 6], &bonds)
            .atoms
            .iter()
            .all(|v| !v)
    );
}

#[test]
fn fused_naphthalene_circuits_are_aromatic() {
    use BondOrder::{Double as D, Single as S};
    let bonds = [
        (0, 1, D),
        (1, 2, S),
        (2, 3, D),
        (3, 4, S),
        (4, 5, D),
        (5, 0, S),
        (4, 6, S),
        (6, 7, D),
        (7, 8, S),
        (8, 9, D),
        (9, 5, S),
    ];
    assert!(
        perceive(&[Element::CARBON; 10], &[Some(0); 10], &bonds)
            .atoms
            .iter()
            .all(|v| *v)
    );
}

#[test]
fn charged_carbon_supplies_two_or_zero_pi_electrons() {
    for (orders, charge) in [
        (
            vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
            -1,
        ),
        (
            vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Single,
            ],
            1,
        ),
    ] {
        let mut charges = vec![Some(0); orders.len()];
        charges[0] = Some(charge);
        let result = ring(&vec![Element::CARBON; orders.len()], &orders, &charges);
        assert!(result.atoms.iter().all(|value| *value));
    }
}

#[test]
fn acyclic_conjugation_does_not_make_aromatic_atoms() {
    let bonds = [
        (0, 1, BondOrder::Double),
        (1, 2, BondOrder::Single),
        (2, 3, BondOrder::Double),
    ];
    assert!(
        perceive(&[Element::CARBON; 4], &[Some(0); 4], &bonds)
            .atoms
            .iter()
            .all(|value| !value)
    );
}
