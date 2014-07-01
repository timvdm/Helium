import helium
import unittest
import os

class TestBinaryInputFile(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME', 'c:/') + '/Helium/data')

    def test_depict(self):
        SMILES = helium.Smiles()
        mol = helium.Molecule()
        self.assertTrue(SMILES.read('c1ccccc1', mol))

        coords = helium.generate_diagram(mol)

        painter = helium.SVGPainter()
        depict = helium.Depict(painter)

        rings = helium.relevant_cycles(mol)
        depict.drawMolecule(mol, rings, coords)

        print painter.output()



if __name__ == '__main__':
    unittest.main()
