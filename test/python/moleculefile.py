import helium
import unittest
import os

class TestBinaryInputFile(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def test_molecule_file(self):
        f = helium.MoleculeFile()
        self.assertFalse(f.load('foo'))

        f.load(self.datadir + '/1K.hel')
        self.assertEqual(1000, f.numMolecules())

        mol = helium.Molecule()
        self.assertTrue(f.readMolecule(mol))
        self.assertEqual(8, mol.numAtoms())
        self.assertTrue(f.readMolecule(998, mol))
        self.assertEqual(29, mol.numAtoms())

        f.close()

        f = helium.MoleculeFile(self.datadir + '/1K.hel')
        self.assertTrue(f.readMolecule(0, mol))
        self.assertEqual(8, mol.numAtoms())

    def test_memory_mapped_molecule_file(self):
        f = helium.MemoryMappedMoleculeFile()
        self.assertFalse(f.load('foo'))

        f.load(self.datadir + '/1K.hel')
        self.assertEqual(1000, f.numMolecules())

        mol = helium.Molecule()
        self.assertTrue(f.readMolecule(0, mol))
        self.assertEqual(8, mol.numAtoms())

        f = helium.MemoryMappedMoleculeFile(self.datadir + '/1K.hel')
        self.assertEqual(1000, f.numMolecules())


if __name__ == '__main__':
    unittest.main()
