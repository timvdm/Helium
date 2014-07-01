import helium
import unittest
import os

SMILES = helium.Smiles()

class BFSTestVisitor(helium.BFSVisitor):

    output = ''

    def initialize(self, mol):
        self.output += 'initialize()\n'

    def component(self, i):
        self.output += 'component(' + str(i) + ')\n'

    def depth(self, d):
        self.output += 'depth(' + str(d) + ')\n'

    def atom(self, mol, prev, atom):
        self.output += 'atom(' + str(atom.index()) + ')\n'

    def bond(self, mol, prev, bond):
        self.output += 'bond(' + str(bond.index()) + ')\n'

    def backtrack(self, mol, atom):
        self.output += 'backtrack(' + str(atom.index()) + ')\n'

    def back_bond(self, mol, bond):
        self.output += 'back_bond(' + str(bond.index()) + ')\n'


class TestBFS(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME', 'c:/') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def compare_file(self, filename, output):
        f = open(self.datadir + '/' + filename)
        lines = f.read()
        self.assertEqual(output, lines)

    def test_bfs_1(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        visitor = helium.BFSDebugVisitor()
        helium.breadth_first_search(mol, visitor)
        self.compare_file('bfs1.log', visitor.output)

    def test_bfs_2(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[4] = False
        atomMask[5] = False

        visitor = helium.BFSDebugVisitor()
        helium.breadth_first_search_mask(mol, visitor, atomMask)
        self.compare_file('bfs2.log', visitor.output)

    def test_bfs_3(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[4] = False
        atomMask[5] = False
        bondMask = [True] * mol.numBonds()
        bondMask[3] = False
        bondMask[4] = False
        bondMask[5] = False

        visitor = helium.BFSDebugVisitor()
        helium.breadth_first_search_mask(mol, visitor, atomMask, bondMask)
        self.compare_file('bfs3.log', visitor.output)

    def test_bfs_4(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        visitor = helium.BFSDebugVisitor()
        helium.breadth_first_search(mol, mol.atom(5), visitor)
        self.compare_file('bfs4.log', visitor.output)

    def test_bfs_5(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[0] = False

        visitor = helium.BFSDebugVisitor()
        helium.breadth_first_search_mask(mol, mol.atom(5), visitor, atomMask)
        self.compare_file('bfs5.log', visitor.output)

    def test_bfs_6(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[0] = False
        bondMask = [True] * mol.numBonds()
        bondMask[0] = False
        bondMask[6] = False
        bondMask[7] = False

        visitor = helium.BFSDebugVisitor()
        helium.breadth_first_search_mask(mol, mol.atom(5), visitor, atomMask, bondMask)
        self.compare_file('bfs6.log', visitor.output)

    def test_atom_order_visitor(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        visitor = helium.BFSAtomOrderVisitor()
        helium.breadth_first_search(mol, visitor)
        self.assertListEqual([0, 1, 2], visitor.atoms)

    def test_bond_order_visitor(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        visitor = helium.BFSBondOrderVisitor()
        helium.breadth_first_search(mol, visitor)
        self.assertListEqual([0, 1], visitor.bonds)

    def test_closure_recorder_visitor(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1', mol)

        visitor = helium.BFSClosureRecorderVisitor()
        helium.breadth_first_search(mol, visitor)
        self.assertListEqual([3], visitor.back_bonds)

    def test_custom_visitor(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        visitor = BFSTestVisitor()
        helium.breadth_first_search(mol, visitor)
        self.compare_file('bfs1.log', visitor.output)


if __name__ == '__main__':
    unittest.main()
