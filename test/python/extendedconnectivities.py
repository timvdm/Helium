import helium
import unittest

SMILES = helium.Smiles()

class TestAtomInvariant(helium.AtomInvariant):

    def __call__(self, mol, atom):
        return atom.element()

class TestExtendedConnectivities(unittest.TestCase):

    def test_extended_connectivities(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        ec = helium.extended_connectivities(mol, helium.DefaultAtomInvariant())
        self.assertTrue(isinstance(ec, list))
        self.assertEqual(3, len(ec))

        ec = helium.extended_connectivities(mol, TestAtomInvariant())
        self.assertTrue(isinstance(ec, list))
        self.assertEqual(3, len(ec))


if __name__ == '__main__':
    unittest.main()
