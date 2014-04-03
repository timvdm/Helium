import helium
import unittest

SMILES = helium.Smiles()

class TestAtomInvariant(helium.AtomInvariant):

    def __call__(self, mol, atom):
        return atom.element()

class TestBondInvariant(helium.BondInvariant):

    def __call__(self, mol, bond):
        return bond.order()

class TestCanonical(unittest.TestCase):

    def test_canonical(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        atomInvariant1 = helium.DefaultAtomInvariant()
        atomInvariant2 = TestAtomInvariant()
        bondInvariant1 = helium.DefaultBondInvariant()
        bondInvariant2 = TestBondInvariant()

        ec = helium.extended_connectivities(mol, atomInvariant1)

        can = helium.canonicalize_component(mol, ec, atomInvariant1, bondInvariant1)
        self.assertTrue(isinstance(can, tuple))
        self.assertTrue(isinstance(can[0], list))
        self.assertTrue(isinstance(can[1], list))
        self.assertEqual(3, len(can[0]))

        can = helium.canonicalize_component(mol, ec, atomInvariant2, bondInvariant1)
        self.assertTrue(isinstance(can, tuple))
        self.assertTrue(isinstance(can[0], list))
        self.assertTrue(isinstance(can[1], list))
        self.assertEqual(3, len(can[0]))

        can = helium.canonicalize_component(mol, ec, atomInvariant1, bondInvariant2)
        self.assertTrue(isinstance(can, tuple))
        self.assertTrue(isinstance(can[0], list))
        self.assertTrue(isinstance(can[1], list))
        self.assertEqual(3, len(can[0]))

        can = helium.canonicalize_component(mol, ec, atomInvariant2, bondInvariant2)
        self.assertTrue(isinstance(can, tuple))
        self.assertTrue(isinstance(can[0], list))
        self.assertTrue(isinstance(can[1], list))
        self.assertEqual(3, len(can[0]))


if __name__ == '__main__':
    unittest.main()
