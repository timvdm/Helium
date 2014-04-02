import helium
import unittest

class TestFingerprints(unittest.TestCase):

    def test_path_fingerprint(self):
        mol = helium.Molecule()
        fp = helium.path_fingerprint(mol)
        fp = helium.path_fingerprint(mol, 7)
        fp = helium.path_fingerprint(mol, 7, 16)
        fp = helium.path_fingerprint(mol, 7, 16, 1021)
        self.assertEqual(16, fp.numWords)

    def test_tree_fingerprint(self):
        mol = helium.Molecule()
        fp = helium.tree_fingerprint(mol)
        fp = helium.tree_fingerprint(mol, 7)
        fp = helium.tree_fingerprint(mol, 7, 16)
        fp = helium.tree_fingerprint(mol, 7, 16, 1021)
        self.assertEqual(16, fp.numWords)

    def test_subgraph_fingerprint(self):
        mol = helium.Molecule()
        fp = helium.subgraph_fingerprint(mol)
        fp = helium.subgraph_fingerprint(mol, 7)
        fp = helium.subgraph_fingerprint(mol, 7, 16)
        fp = helium.subgraph_fingerprint(mol, 7, 16, 1021)
        self.assertEqual(16, fp.numWords)

if __name__ == '__main__':
    unittest.main()
