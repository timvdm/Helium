import helium
import unittest

SMILES = helium.Smiles()

class TestSubgraphCallback(helium.SubgraphCallback):

    subgraphs = []

    def __call__(self, subgraph):
        self.subgraphs.append((subgraph.atoms, subgraph.bonds))

class TestEnumerateSubgraphs(unittest.TestCase):

    def test_subgraph(self):
        subgraph = helium.Subgraph(3, 2)
        self.assertTrue(isinstance(subgraph.hashable(), list))
        self.assertTrue(isinstance(subgraph.atoms, list))
        self.assertTrue(isinstance(subgraph.bonds, list))

    def test_enumerate_paths(self):
        mol = helium.Molecule()
        SMILES.read('CC(C)C', mol)

        callback = TestSubgraphCallback()

        helium.enumerate_subgraphs(mol, callback, 7)
        self.assertTrue(isinstance(callback.subgraphs, list))
        self.assertEqual(11, len(callback.subgraphs))

        callback.subgraphs = []
        helium.enumerate_subgraphs(mol, callback, 3)
        self.assertTrue(isinstance(callback.subgraphs, list))
        self.assertEqual(10, len(callback.subgraphs))


if __name__ == '__main__':
    unittest.main()
