import helium
import unittest

SMILES = helium.Smiles()

class TestSmiles(unittest.TestCase):

    def test_read_valid(self):
        mol = helium.Molecule()
        self.assertTrue(SMILES.read('CCC', mol))
        self.assertFalse(SMILES.error())

    def test_read_invalid(self):
        mol = helium.Molecule()
        self.assertFalse(SMILES.read('dfgd', mol))
        self.assertTrue(SMILES.error())
        self.assertNotEqual(0, len(str(SMILES.error())))

    def test_write(self):
        mol = helium.Molecule()
        self.assertTrue(SMILES.read('C=C', mol))
        self.assertEqual('C=C', SMILES.write(mol))
        self.assertEqual('CC', SMILES.write(mol, helium.Smiles.Flags.None))

    def test_canonical(self):
        mol1 = helium.Molecule()
        mol2 = helium.Molecule()
        self.assertTrue(SMILES.read('C=O', mol1))
        self.assertTrue(SMILES.read('O=C', mol2))

        self.assertEqual(SMILES.writeCanonical(mol1), SMILES.writeCanonical(mol2))
        self.assertEqual(SMILES.write(mol1, [0, 1]), SMILES.write(mol2, [1, 0]))
        self.assertEqual('CO', SMILES.write(mol1, [0, 1], helium.Smiles.Flags.None))
        self.assertEqual('CO', SMILES.write(mol2, [1, 0], helium.Smiles.Flags.None))


if __name__ == '__main__':
    unittest.main()
