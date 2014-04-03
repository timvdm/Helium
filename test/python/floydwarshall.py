import helium
import unittest

class TestFloydWarshall(unittest.TestCase):

    def test_floyd_warshall(self):
        mol = helium.Molecule()

        a1 = mol.addAtom()
        a2 = mol.addAtom()
        a3 = mol.addAtom()
        a4 = mol.addAtom()
        a5 = mol.addAtom()
        a6 = mol.addAtom()
        a7 = mol.addAtom()

        mol.addBond(a1, a2)
        mol.addBond(a2, a3)
        mol.addBond(a3, a4)
        mol.addBond(a4, a5)
        mol.addBond(a5, a6)
        mol.addBond(a6, a7)
        mol.addBond(a7, a1)
        mol.addBond(a2, a6)

        D = helium.floyd_warshall(mol)
        self.assertTrue(isinstance(D, helium.DistanceMatrix))
        self.assertEqual(7, D.dim())
        self.assertEqual(0, D(0, 0))
        self.assertEqual(1, D(0, 1))
        self.assertEqual(2, D(0, 2))
        self.assertEqual(3, D(0, 3))
        self.assertEqual(4, D(0, 4))
        self.assertEqual(2, D(0, 5))
        self.assertEqual(3, D(0, 6))

        self.assertTrue(isinstance(D.infinity(), int))


if __name__ == '__main__':
    unittest.main()
