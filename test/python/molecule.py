import helium
import unittest

SMILES = helium.Smiles()

class TestMolecule(unittest.TestCase):

    def test_molecule(self):
        mol = helium.Molecule()
        SMILES.read('CC(C)N', mol)

        self.assertEqual(4, mol.numAtoms())
        self.assertEqual(3, mol.numBonds())

        atoms = []
        for atom in mol.atoms():
            atoms.append(atom)
        self.assertListEqual([mol.atom(0), mol.atom(1), mol.atom(2), mol.atom(3)], atoms)

        bonds = []
        for bond in mol.bonds():
            bonds.append(bond)
        self.assertListEqual([mol.bond(0), mol.bond(1), mol.bond(2)], bonds)

    def test_atom(self):
        mol = helium.Molecule()
        SMILES.read('CC(C)N', mol)

        atom1 = mol.atom(0)
        atom2 = mol.atom(1)
        atom3 = mol.atom(2)
        atom4 = mol.atom(3)

        self.assertEqual(0, atom1.index())
        self.assertEqual(1, atom2.index())
        self.assertEqual(2, atom3.index())
        self.assertEqual(3, atom4.index())

        self.assertEqual(6, atom1.element())
        self.assertEqual(6, atom2.element())
        self.assertEqual(6, atom3.element())
        self.assertEqual(7, atom4.element())

        self.assertEqual(1, atom1.degree())
        self.assertEqual(3, atom2.degree())
        self.assertEqual(1, atom3.degree())
        self.assertEqual(1, atom4.degree())

        indices = []
        for nbr in atom2.nbrs():
            indices.append(nbr.index())
        self.assertListEqual([0, 2, 3], indices)

        indices = []
        for bond in atom2.bonds():
            indices.append(bond.index())
        self.assertListEqual([0, 1, 2], indices)

    def test_bond(self):
        mol = helium.Molecule()
        SMILES.read('CC(C)=O', mol)

        atom1 = mol.atom(0)
        atom2 = mol.atom(1)
        atom3 = mol.atom(2)
        atom4 = mol.atom(3)

        bond1 = mol.bond(0)
        bond2 = mol.bond(1)
        bond3 = mol.bond(2)

        self.assertEqual(0, bond1.index())
        self.assertEqual(1, bond2.index())
        self.assertEqual(2, bond3.index())

        self.assertEqual(1, bond1.order())
        self.assertEqual(1, bond2.order())
        self.assertEqual(2, bond3.order())

        self.assertEqual(atom1, bond1.source())
        self.assertEqual(atom2, bond2.source())
        self.assertEqual(atom2, bond3.source())

        self.assertEqual(atom2, bond1.target())
        self.assertEqual(atom3, bond2.target())
        self.assertEqual(atom4, bond3.target())

        self.assertEqual(atom1, bond1.other(atom2))
        self.assertEqual(atom2, bond2.other(atom3))
        self.assertEqual(atom2, bond3.other(atom4))

    def test_atom_editing(self):
        mol = helium.Molecule()
        atom = mol.addAtom()

        self.assertEqual(0, atom.index())

        atom.setElement(6)
        self.assertEqual(6, atom.element())
        atom.setElement(7)
        self.assertEqual(7, atom.element())

        atom.setAromatic(True)
        self.assertTrue(atom.isAromatic())
        atom.setAromatic(False)
        self.assertFalse(atom.isAromatic())

        atom.setMass(12)
        self.assertEqual(12, atom.mass())
        atom.setMass(13)
        self.assertEqual(13, atom.mass())

        atom.setHydrogens(4)
        self.assertEqual(4, atom.hydrogens())
        atom.setHydrogens(3)
        self.assertEqual(3, atom.hydrogens())

        atom.setCharge(1)
        self.assertEqual(1, atom.charge())
        atom.setCharge(-1)
        self.assertEqual(-1, atom.charge())

    def test_bond_editing(self):
        mol = helium.Molecule()
        source = mol.addAtom()
        target = mol.addAtom()
        bond = mol.addBond(source, target)

        bond.setAromatic(True)
        self.assertTrue(bond.isAromatic())
        bond.setAromatic(False)
        self.assertFalse(bond.isAromatic())

        bond.setOrder(1)
        self.assertEqual(1, bond.order())
        bond.setOrder(2)
        self.assertEqual(2, bond.order())

    def test_molecule_editing(self):
        mol = helium.Molecule()

        atom1 = mol.addAtom()
        atom2 = mol.addAtom()
        atom3 = mol.addAtom()
        atom4 = mol.addAtom()

        bond1 = mol.addBond(atom1, atom2)
        bond2 = mol.addBond(atom2, atom3)
        bond3 = mol.addBond(atom3, atom4)

        self.assertEqual(4, mol.numAtoms())
        self.assertEqual(3, mol.numBonds())

        mol.removeBond(bond3)

        self.assertEqual(4, mol.numAtoms())
        self.assertEqual(2, mol.numBonds())

        mol.removeAtom(atom1)

        self.assertEqual(3, mol.numAtoms())
        self.assertEqual(1, mol.numBonds())

        mol.clear()

        self.assertEqual(0, mol.numAtoms())
        self.assertEqual(0, mol.numBonds())

    def test_copy_ctor(self):
        mol1 = helium.Molecule()

        SMILES.read('CC', mol1)

        mol2 = helium.Molecule(mol1)

        self.assertEqual(2, mol1.numAtoms())
        self.assertEqual(2, mol2.numAtoms())

        mol2.clear()

        self.assertEqual(2, mol1.numAtoms())
        self.assertEqual(0, mol2.numAtoms())

    def test_get_bond(self):
        mol = helium.Molecule()
        SMILES.read('C=C', mol)

        bond = mol.bond(0)
        bond = mol.bond(mol.atom(0), mol.atom(1))

if __name__ == '__main__':
    unittest.main()
