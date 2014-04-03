import helium
import unittest

SMILES = helium.Smiles()

class TestEnumeratePaths(unittest.TestCase):

    def test_enumerate_paths(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        paths = helium.enumerate_paths(mol, 7)
        self.assertTrue(isinstance(paths, list))
        self.assertEqual(6, len(paths))
        self.assertListEqual([0], paths[0])
        self.assertListEqual([1], paths[1])
        self.assertListEqual([2], paths[2])
        self.assertListEqual([0, 1], paths[3])
        self.assertListEqual([1, 2], paths[4])
        self.assertListEqual([0, 1, 2], paths[5])

        paths = helium.enumerate_paths(mol, 2)
        self.assertEqual(5, len(paths))
        self.assertListEqual([0], paths[0])
        self.assertListEqual([1], paths[1])
        self.assertListEqual([2], paths[2])
        self.assertListEqual([0, 1], paths[3])
        self.assertListEqual([1, 2], paths[4])


if __name__ == '__main__':
    unittest.main()
