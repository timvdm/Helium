import helium
import unittest

SMILES = helium.Smiles()

class TestRings(unittest.TestCase):

    def test_cyclomatic_number(self):
        mol = helium.Molecule()
        SMILES.read('CC1CCCCC1C', mol)

        self.assertEqual(1, helium.cyclomatic_number(mol, 1))
        self.assertEqual(1, helium.cyclomatic_number(mol))

    def test_cycle_membership(self):
        mol = helium.Molecule()
        SMILES.read('CC1CC1C', mol)

        atoms, bonds = helium.cycle_membership(mol)

        self.assertListEqual([False, True, True, True, False], atoms)
        self.assertListEqual([False, True, True, True, False], bonds)

    def test_relevant_cycles(self):
        mol = helium.Molecule()
        SMILES.read('c1ccccc1-c1ccccc1', mol)

        rings = helium.relevant_cycles(mol)
        self.assertEqual(2, len(rings))

        atoms, bonds = helium.cycle_membership(mol)
        rings = helium.relevant_cycles(mol, 1, atoms, bonds)
        self.assertEqual(2, len(rings))

    def test_ringset(self):
        mol = helium.Molecule()
        SMILES.read('CCC1CC1', mol)
        rings = helium.relevant_cycles(mol)
        self.assertTrue(isinstance(rings, helium.RingSet))

        self.assertEqual(1, rings.size())
        self.assertEqual(1, len(rings))

        self.assertEqual(1, len(rings.rings()))
        self.assertTrue(isinstance(rings.ring(0), helium.Ring))

        self.assertFalse(rings.isAtomInRing(mol.atom(0)))
        self.assertTrue(rings.isAtomInRing(mol.atom(2)))

        self.assertFalse(rings.isBondInRing(mol.bond(0)))
        self.assertTrue(rings.isBondInRing(mol.bond(2)))

        self.assertFalse(rings.isAtomInRingSize(mol.atom(2), 4))
        self.assertTrue(rings.isAtomInRingSize(mol.atom(2), 3))

        self.assertFalse(rings.isBondInRingSize(mol.bond(2), 4))
        self.assertTrue(rings.isBondInRingSize(mol.bond(2), 3))

        self.assertEqual(2, rings.numRingBonds(mol.atom(2)))
        self.assertEqual(2, rings.numRingNbrs(mol.atom(2)))
        self.assertEqual(1, rings.numRings(mol.atom(2)))

    def test_ring(self):
        mol = helium.Molecule()
        SMILES.read('CCC1CC1', mol)
        rings = helium.relevant_cycles(mol)
        ring = rings.ring(0)

        self.assertEqual(3, ring.size())
        self.assertEqual(3, len(ring))
        self.assertEqual(3, len(ring.atoms()))
        self.assertEqual(3, len(ring.bonds()))

        self.assertEqual(4, ring.atom(0).index())
        self.assertEqual(3, ring.atom(1).index())
        self.assertEqual(2, ring.atom(2).index())

        self.assertEqual(3, ring.bond(0).index())
        self.assertEqual(2, ring.bond(1).index())
        self.assertEqual(4, ring.bond(2).index())

        self.assertFalse(ring.containsAtom(mol.atom(0)))
        self.assertTrue(ring.containsAtom(mol.atom(2)))

        self.assertFalse(ring.containsBond(mol.bond(0)))
        self.assertTrue(ring.containsBond(mol.bond(2)))

if __name__ == '__main__':
    unittest.main()
