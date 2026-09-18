use super::*;

#[test]
fn optimization_reads_all_atoms_energy_and_bohr_units() {
    let source = " INPUT CARD> $CONTRL RUNTYP=OPTIMIZE $END\n\
 TOTAL NUMBER OF ATOMS = 2\n\
1NSERCH= 0\n COORDINATES OF ALL ATOMS ARE (BOHR)\n ATOM CHARGE X Y Z\n -----\n\
C 6.0 0 0 0\nH 1.0 2 0 0\n NSERCH= 0 ENERGY= -40.25\n\
1NSERCH= 1\n COORDINATES OF ALL ATOMS ARE (ANGS)\n ATOM CHARGE X Y Z\n -----\n\
C 6.0 0 0 0\nH 1.0 1.1 0 0\n NSERCH= 1 ENERGY= -40.50\n";
    let trajectory = parse_gamess_output(source)
        .unwrap_or_else(|error| panic!("optimization parse failed: {error}"));
    assert_eq!(trajectory.frames.len(), 2);
    assert!((trajectory.frames[0].positions[1][0] - 1.058_354_4).abs() < 1.0e-6);
    assert_eq!(trajectory.frames[1].energy, Some(-40.5));
    assert_eq!(trajectory.atoms[0].nuclear_charge, Some(6.0));
}

#[test]
fn surface_scan_retains_coordinates_energy_and_stable_names() {
    let source = "$CONTRL RUNTYP=SURFACE $END\nTOTAL NUMBER OF ATOMS = 2\n\
 COORD 1= 0.0 COORD 2= -0.1\n HAS ENERGY VALUE -75.0\n\
O 0 0 0\nH 1 0 0\n\
 COORD 1= 0.2 COORD 2= -0.1\n HAS ENERGY VALUE -75.2\n\
O 0 0 0\nH 1.2 0 0\n";
    let trajectory =
        parse_gamess_output(source).unwrap_or_else(|error| panic!("surface parse failed: {error}"));
    assert_eq!(trajectory.run_type, GamessRunType::Surface);
    assert_eq!(trajectory.frames[1].surface_coordinates, Some([0.2, -0.1]));
    let timestep = trajectory.frames[1].to_timestep(1);
    assert_eq!(
        timestep.data.get("gamess.energy_hartree"),
        Some(&FrameValue::Float(-75.2))
    );
}

#[test]
fn changed_atom_name_is_rejected_as_topology_drift() {
    let source = "$CONTRL RUNTYP=SURFACE $END\nTOTAL NUMBER OF ATOMS = 1\n\
COORD 1= 0 COORD 2= 0\nHAS ENERGY VALUE -1\nH 0 0 0\n\
COORD 1= 1 COORD 2= 0\nHAS ENERGY VALUE -2\nO 0 0 0\n";
    assert_eq!(parse_gamess_output(source), Err(GamessError::TopologyDrift));
}
