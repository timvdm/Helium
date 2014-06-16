import helium
import unittest

SMILES = helium.Smiles()

class TestInvariants(unittest.TestCase):

    def test_default_atom_invariant(self):
        mol = helium.Molecule()
        SMILES.read('CNO', mol)

        invariant = helium.DefaultAtomInvariant()
        self.assertEqual(2110012006, invariant(mol, mol.atom(0)))

        invariant = helium.DefaultAtomInvariant(helium.DefaultAtomInvariant.Invariants.Element)
        self.assertEqual(6, invariant(mol, mol.atom(0)))

        invariant = helium.DefaultAtomInvariant(helium.DefaultAtomInvariant.Invariants.Mass)
        self.assertEqual(12000, invariant(mol, mol.atom(0)))

        ec = helium.extended_connectivities(mol, invariant)
        self.assertTrue(isinstance(ec, list))
        self.assertEqual(3, len(ec))

    def test_default_bond_invariant(self):
        mol = helium.Molecule()
        SMILES.read('CN=O:C', mol)

        invariant = helium.DefaultBondInvariant()
        self.assertEqual(1, invariant(mol, mol.bond(0)))

        invariant = helium.DefaultBondInvariant(helium.DefaultBondInvariant.Invariants.Order)
        self.assertEqual(2, invariant(mol, mol.bond(1)))

        invariant = helium.DefaultBondInvariant(helium.DefaultBondInvariant.Invariants.Aromatic)
        self.assertEqual(100, invariant(mol, mol.bond(2)))


if __name__ == '__main__':
    unittest.main()
