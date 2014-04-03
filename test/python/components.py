import helium
import unittest

SMILES = helium.Smiles()

class TestComponents(unittest.TestCase):

    def test_components(self):
        mol = helium.Molecule()

        SMILES.read('CC(C)C', mol)
        components = helium.connected_bond_components(mol)
        self.assertTrue(isinstance(components, list))
        self.assertEqual(3, len(components))
        self.assertEqual(1, len(set(components)))

        SMILES.read('CC.CC', mol)
        components = helium.connected_bond_components(mol)
        self.assertTrue(isinstance(components, list))
        self.assertEqual(2, len(components))
        self.assertEqual(2, len(set(components)))
        self.assertEqual(0, components[0])
        self.assertEqual(1, components[1])

        SMILES.read('C1CC1.CC', mol)
        components = helium.connected_bond_components(mol)
        self.assertTrue(isinstance(components, list))
        self.assertEqual(4, len(components))
        self.assertEqual(2, len(set(components)))
        self.assertEqual(0, components[0])
        self.assertEqual(0, components[1])
        self.assertEqual(0, components[2])
        self.assertEqual(1, components[3])

        SMILES.read('CC(C)C', mol)
        components = helium.connected_atom_components(mol)
        self.assertTrue(isinstance(components, list))
        self.assertEqual(4, len(components))
        self.assertEqual(1, len(set(components)))

        SMILES.read('CC.CC', mol)
        components = helium.connected_atom_components(mol)
        self.assertTrue(isinstance(components, list))
        self.assertEqual(4, len(components))
        self.assertEqual(2, len(set(components)))
        self.assertEqual(0, components[0])
        self.assertEqual(0, components[1])
        self.assertEqual(1, components[2])
        self.assertEqual(1, components[3])

        SMILES.read('C1CC1.CC', mol)
        components = helium.connected_atom_components(mol)
        self.assertTrue(isinstance(components, list))
        self.assertEqual(5, len(components))
        self.assertEqual(2, len(set(components)))
        self.assertEqual(0, components[0])
        self.assertEqual(0, components[1])
        self.assertEqual(0, components[2])
        self.assertEqual(1, components[3])
        self.assertEqual(1, components[4])

        SMILES.read('CC(C)C', mol)
        self.assertEqual(1, helium.num_connected_components(mol))

        SMILES.read('CC.CC', mol)
        self.assertEqual(2, helium.num_connected_components(mol))

        SMILES.read('C1CC1.CC', mol)
        self.assertEqual(2, helium.num_connected_components(mol))

if __name__ == '__main__':
    unittest.main()
