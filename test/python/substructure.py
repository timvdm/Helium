import helium
import unittest

SMILES = helium.Smiles()

class TestSubstructure(unittest.TestCase):

    def test_substructure(self):
        mol = helium.Molecule()
        SMILES.read('Oc1cc(CC)ccc1N', mol)

        atoms = [True] * mol.numAtoms()
        atoms[0] = False
        atoms[4] = False
        atoms[5] = False
        atoms[9] = False

        bonds = [True] * mol.numBonds()
        bonds[0] = False
        bonds[3] = False
        bonds[4] = False
        bonds[9] = False

        sub = helium.Substructure(mol, atoms, bonds)

        self.assertEqual('c1ccccc1', SMILES.write(sub))

    def test_exceptions(self):
        mol = helium.Molecule()
        SMILES.read('CCCCC', mol)

        # test invalid atoms
        with self.assertRaises(RuntimeError):
            helium.Substructure(mol, [], [False] * mol.numBonds())

        # test invalid bonds
        with self.assertRaises(RuntimeError):
            helium.Substructure(mol, [False] * mol.numAtoms(), [])




if __name__ == '__main__':
    unittest.main()
