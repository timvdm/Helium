import helium
import unittest

SMILES = helium.Smiles()

class TestSubstructure(unittest.TestCase):

    def test_substructure(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1C', mol)

        sub = helium.Substructure(mol)



if __name__ == '__main__':
    unittest.main()
