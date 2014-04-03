import helium
import unittest

SMILES = helium.Smiles()

class TestDijkstra(unittest.TestCase):

    def test_dijkstra(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCC2C1CCC2', mol)

        d = helium.Dijkstra(mol, mol.atom(0))

        self.assertEqual(0, d.distance(mol.atom(0)))
        self.assertEqual(1, d.distance(mol.atom(1)))
        self.assertEqual(2, d.distance(mol.atom(2)))
        self.assertEqual(3, d.distance(mol.atom(3)))
        self.assertEqual(2, d.distance(mol.atom(4)))
        self.assertEqual(1, d.distance(mol.atom(5)))
        self.assertEqual(2, d.distance(mol.atom(6)))
        self.assertEqual(3, d.distance(mol.atom(7)))
        self.assertEqual(3, d.distance(mol.atom(8)))

        path = d.path(mol.atom(0))
        self.assertEqual(1, len(path))
        self.assertEqual(0, path[0].index())

        path = d.path(mol.atom(1))
        self.assertEqual(2, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(1, path[1].index())

        path = d.path(mol.atom(2))
        self.assertEqual(3, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(1, path[1].index())
        self.assertEqual(2, path[2].index())

        path = d.path(mol.atom(3))
        self.assertEqual(4, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(1, path[1].index())
        self.assertEqual(2, path[2].index())
        self.assertEqual(3, path[3].index())

        path = d.path(mol.atom(4))
        self.assertEqual(3, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(5, path[1].index())
        self.assertEqual(4, path[2].index())

        path = d.path(mol.atom(5))
        self.assertEqual(2, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(5, path[1].index())

        path = d.path(mol.atom(6))
        self.assertEqual(3, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(5, path[1].index())
        self.assertEqual(6, path[2].index())

        path = d.path(mol.atom(7))
        self.assertEqual(4, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(5, path[1].index())
        self.assertEqual(6, path[2].index())
        self.assertEqual(7, path[3].index())

        path = d.path(mol.atom(8))
        self.assertEqual(4, len(path))
        self.assertEqual(0, path[0].index())
        self.assertEqual(5, path[1].index())
        self.assertEqual(4, path[2].index())
        self.assertEqual(8, path[3].index())

        self.assertTrue(isinstance(d.infinity(), int))

if __name__ == '__main__':
    unittest.main()
