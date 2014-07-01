import helium
import unittest
import os

SMILES = helium.Smiles()

class TestSimilarity(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME', 'c:/') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def test_brute_force_similarity_search(self):
        storage = helium.InMemoryRowMajorFingerprintStorage()
        storage.load(self.datadir + '/1K_fp_path_7_1024_row.hel')

        mol = helium.Molecule()
        self.assertTrue(SMILES.read('c1ccccc1', mol))
        query = helium.path_fingerprint(mol)

        similar = helium.brute_force_similarity_search(query, storage, 0.8)
        self.assertTrue(isinstance(similar, list))
        self.assertEqual(1, len(similar))
        self.assertTrue(isinstance(similar[0], tuple))
        i, tanimoto = similar[0]
        self.assertEqual(227, i)
        self.assertTrue(tanimoto > 0.8)

    def test_similarity_search_index(self):
        storage = helium.InMemoryRowMajorFingerprintStorage()
        storage.load(self.datadir + '/1K_fp_path_7_1024_row.hel')

        mol = helium.Molecule()
        self.assertTrue(SMILES.read('c1ccccc1', mol))
        query = helium.path_fingerprint(mol)

        index = helium.SimilaritySearchIndex(storage, 3)
        similar = index.search(query, 0.8)
        similar = index.search(query, 0.8)

        self.assertTrue(isinstance(similar, list))
        self.assertEqual(1, len(similar))
        self.assertTrue(isinstance(similar[0], tuple))
        i, tanimoto = similar[0]
        self.assertEqual(227, i)
        self.assertTrue(tanimoto > 0.8)

        similar = index.search(query, 0.4, 10)
        self.assertEqual(10, len(similar))


if __name__ == '__main__':
    unittest.main()
