use crate::structure::PyStructure;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use std::ffi::CString;
use std::path::PathBuf;

const MOCKS: &str = r"
import sys
import types

class Element:
    @staticmethod
    def getByAtomicNumber(number):
        return number

class Topology:
    def __init__(self):
        self.chains = []
        self.residues = []
        self.atoms = []
        self.bonds = []
    def addChain(self, id):
        chain = {'id': id}
        self.chains.append(chain)
        return chain
    def addResidue(self, name, chain, **fields):
        residue = {'name': name, 'chain': chain, **fields}
        self.residues.append(residue)
        return residue
    def addAtom(self, name, element, residue, **fields):
        atom = {'name': name, 'element': element, 'residue': residue, **fields}
        self.atoms.append(atom)
        return atom
    def addBond(self, atom_a, atom_b, **fields):
        self.bonds.append((atom_a, atom_b, fields))

openmm = types.ModuleType('openmm')
app = types.ModuleType('openmm.app')
app.Topology = Topology
app.element = types.SimpleNamespace(Element=Element)
unit = types.ModuleType('openmm.unit')
unit.angstrom = 1.0
openmm.app = app
openmm.unit = unit
sys.modules['openmm'] = openmm
sys.modules['openmm.app'] = app
sys.modules['openmm.unit'] = unit

class Atom:
    def __init__(self, **fields):
        self.__dict__.update(fields)

class Bond:
    def __init__(self, atom_a, atom_b, order=None):
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.order = order

class Structure:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.residue_arguments = []
        self.coordinates = None
    def add_atom(self, atom, residue_name, residue_number, chain, insertion_code):
        self.atoms.append(atom)
        self.residue_arguments.append(
            (residue_name, residue_number, chain, insertion_code)
        )

parmed = types.ModuleType('parmed')
parmed.Atom = Atom
parmed.Bond = Bond
parmed.Structure = Structure
sys.modules['parmed'] = parmed

class AseAtoms:
    def __init__(self, **fields):
        self.__dict__.update(fields)

ase = types.ModuleType('ase')
ase.Atoms = AseAtoms
sys.modules['ase'] = ase

class AtomArray:
    def __init__(self, size):
        self.size = size
        self.bonds = None

class BondList:
    def __init__(self, atom_count, bonds):
        self.atom_count = atom_count
        self.values = bonds

biotite = types.ModuleType('biotite')
biotite_structure = types.ModuleType('biotite.structure')
biotite_structure.AtomArray = AtomArray
biotite_structure.BondList = BondList
biotite.structure = biotite_structure
sys.modules['biotite'] = biotite
sys.modules['biotite.structure'] = biotite_structure

class RdAtom:
    def __init__(self, atomic_number):
        self.atomic_number = atomic_number
        self.formal_charge = None
    def SetFormalCharge(self, charge):
        self.formal_charge = charge

class RWMol:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.conformers = []
    def AddAtom(self, atom):
        self.atoms.append(atom)
    def AddBond(self, atom_a, atom_b, order):
        self.bonds.append((atom_a, atom_b, order))
    def GetMol(self):
        return self
    def AddConformer(self, conformer, assign_id):
        self.conformers.append(conformer)

class Conformer:
    def __init__(self, size):
        self.positions = [None] * size
    def SetAtomPosition(self, index, point):
        self.positions[index] = point

class Point3D:
    def __init__(self, x, y, z):
        self.xyz = (x, y, z)

chem = types.ModuleType('rdkit.Chem')
chem.RWMol = RWMol
chem.Atom = RdAtom
chem.Conformer = Conformer
chem.BondType = types.SimpleNamespace(
    SINGLE=1, DOUBLE=2, TRIPLE=3, QUADRUPLE=4, AROMATIC=5, UNSPECIFIED=0
)
geometry = types.ModuleType('rdkit.Geometry')
geometry.Point3D = Point3D
rdkit = types.ModuleType('rdkit')
rdkit.Chem = chem
rdkit.Geometry = geometry
sys.modules['rdkit'] = rdkit
sys.modules['rdkit.Chem'] = chem
sys.modules['rdkit.Geometry'] = geometry

class OstEditor:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.updated = False
    def InsertChain(self, name):
        return {'name': name}
    def AppendResidue(self, chain, name, number=None):
        return {'chain': chain, 'name': name, 'number': number}
    def InsertAtom(self, residue, name, position, element, occupancy, bfactor, hetero):
        atom = (residue, name, position, element, occupancy, bfactor, hetero)
        self.atoms.append(atom)
        return atom
    def Connect(self, atom_a, atom_b):
        self.bonds.append((atom_a, atom_b))
    def UpdateICS(self):
        self.updated = True

class OstEntity:
    def __init__(self):
        self.editor = OstEditor()
    def EditXCS(self, mode):
        return self.editor

class Vec3:
    def __init__(self, x, y, z):
        self.xyz = (x, y, z)

ost = types.ModuleType('ost')
ost_mol = types.ModuleType('ost.mol')
ost_geom = types.ModuleType('ost.geom')
ost_mol.CreateEntity = OstEntity
ost_mol.EditMode = types.SimpleNamespace(BUFFERED_EDIT=1)
ost_geom.Vec3 = Vec3
ost.mol = ost_mol
ost.geom = ost_geom
sys.modules['ost'] = ost
sys.modules['ost.mol'] = ost_mol
sys.modules['ost.geom'] = ost_geom
";

#[test]
fn object_model_adapters_receive_native_projection() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let structure = match pdbiox::read(fixture) {
        Ok(structure) => structure,
        Err(findings) => panic!("native fixture failed to read: {findings:?}"),
    };
    Python::initialize();
    Python::attach(|py| {
        if PyModule::import(py, "numpy").is_err() {
            return;
        }
        install_mock_modules(py);
        let bound = Bound::new(py, PyStructure::new(structure.clone()))
            .expect("structure should bind to Python");

        let openmm = bound
            .call_method0("to_openmm")
            .expect("OpenMM export should succeed");
        let pair = openmm.cast::<PyTuple>().expect("OpenMM export is a pair");
        let topology = pair.get_item(0).expect("topology is present");
        let positions = pair.get_item(1).expect("positions are present");
        let atom_count: usize = topology
            .getattr("atoms")
            .expect("topology atoms")
            .len()
            .expect("atom list length");
        assert_eq!(atom_count, structure.atom_count() as usize);
        assert_eq!(
            positions
                .getattr("shape")
                .expect("coordinate shape")
                .extract::<(usize, usize)>()
                .expect("two-dimensional shape"),
            (structure.atom_count() as usize, 3)
        );

        let ase = bound
            .call_method0("to_ase")
            .expect("ASE export should succeed");
        assert_eq!(
            ase.getattr("numbers")
                .expect("ASE numbers")
                .len()
                .expect("number count"),
            structure.atom_count() as usize
        );

        let biotite = bound
            .call_method0("to_biotite")
            .expect("Biotite export should succeed");
        assert_eq!(
            biotite
                .getattr("size")
                .expect("Biotite size")
                .extract::<usize>()
                .expect("integer size"),
            structure.atom_count() as usize
        );

        let rdkit = bound
            .call_method0("to_rdkit")
            .expect("RDKit export should succeed");
        assert_eq!(
            rdkit
                .getattr("atoms")
                .expect("RDKit atoms")
                .len()
                .expect("atom count"),
            structure.atom_count() as usize
        );

        let openstructure = bound
            .call_method0("to_openstructure")
            .expect("OpenStructure export should succeed");
        assert_eq!(
            openstructure
                .getattr("editor")
                .expect("OpenStructure editor")
                .getattr("atoms")
                .expect("OpenStructure atoms")
                .len()
                .expect("atom count"),
            structure.atom_count() as usize
        );

        let parmed = bound
            .call_method0("to_parmed")
            .expect("ParmEd export should succeed");
        let parmed_atom_count = parmed
            .getattr("atoms")
            .expect("ParmEd atoms")
            .len()
            .expect("atom list length");
        assert_eq!(parmed_atom_count, structure.atom_count() as usize);
        assert_eq!(
            parmed
                .getattr("coordinates")
                .expect("ParmEd coordinates")
                .getattr("shape")
                .expect("coordinate shape")
                .extract::<(usize, usize)>()
                .expect("two-dimensional shape"),
            (structure.atom_count() as usize, 3)
        );
    });
}

#[test]
fn installed_object_model_adapters_when_configured() {
    if std::env::var_os("PDBIOX_EXTERNAL_ADAPTER_TESTS").is_none() {
        return;
    }
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/basic.pdb");
    let structure = match pdbiox::read(fixture) {
        Ok(structure) => structure,
        Err(findings) => panic!("native fixture failed to read: {findings:?}"),
    };
    Python::initialize();
    Python::attach(|py| {
        let bound =
            Bound::new(py, PyStructure::new(structure)).expect("structure should bind to Python");
        for method in [
            "to_gemmi",
            "to_biopython",
            "to_biotite",
            "to_mdanalysis",
            "to_mdtraj",
        ] {
            bound
                .call_method0(method)
                .unwrap_or_else(|error| panic!("{method} failed: {error}"));
        }
    });
}

fn install_mock_modules(py: Python<'_>) {
    let code = CString::new(MOCKS).expect("mock source has no null bytes");
    py.run(&code, None, None)
        .expect("mock external modules should install");
}
