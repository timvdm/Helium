import helium
import unittest
import os

SMILES = helium.Smiles()

class DFSTestVisitor(helium.DFSVisitor):

    output = ''

    def initialize(self, mol):
        self.output += 'initialize()\n'

    def component(self, i):
        self.output += 'component(' + str(i) + ')\n'

    def atom(self, mol, prev, atom):
        self.output += 'atom(' + str(atom.index()) + ')\n'

    def bond(self, mol, prev, bond):
        self.output += 'bond(' + str(bond.index()) + ')\n'

    def backtrack(self, mol, atom):
        self.output += 'backtrack(' + str(atom.index()) + ')\n'

    def back_bond(self, mol, bond):
        self.output += 'back_bond(' + str(bond.index()) + ')\n'


class TestDFS(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME', 'c:/') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def compare_file(self, filename, output):
        f = open(self.datadir + '/' + filename)
        lines = f.read()
        self.assertEqual(output, lines)

    def test_dfs_1(self):
        mol = helium.Molecule()
        SMILES.read('C(CC)(CC)CC', mol)

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search(mol, visitor)
        self.compare_file('dfs1.log', visitor.output)

    def test_dfs_2(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search(mol, visitor)
        self.compare_file('dfs2.log', visitor.output)

    def test_dfs_3(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        visitor = helium.DFSDebugVisitor()
        helium.exhaustive_depth_first_search(mol, mol.atom(0), visitor)
        self.compare_file('dfs3.log', visitor.output)

    def test_dfs_4(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1CC.CC', mol)

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search(mol, visitor)
        self.compare_file('dfs4.log', visitor.output)

    def test_dfs_5(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1CC.CC', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[6] = False
        atomMask[7] = False

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search_mask(mol, visitor, atomMask)
        self.compare_file('dfs5.log', visitor.output)

    def test_dfs_6(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1CC.CC', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[6] = False
        atomMask[7] = False
        bondMask = [True] * mol.numBonds()
        bondMask[6] = False
        bondMask[7] = False
        bondMask[8] = False

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search_mask(mol, visitor, atomMask, bondMask)
        self.compare_file('dfs6.log', visitor.output)

    def test_dfs_7(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1CC.CC', mol)

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search(mol, mol.atom(0), visitor)
        self.compare_file('dfs7.log', visitor.output)

    def test_dfs_8(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1CC.CC', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[6] = False
        atomMask[7] = False

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search_mask(mol, mol.atom(0), visitor, atomMask)
        self.compare_file('dfs8.log', visitor.output)

    def test_dfs_9(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1CC.CC', mol)

        atomMask = [True] * mol.numAtoms()
        atomMask[6] = False
        atomMask[7] = False
        bondMask = [True] * mol.numBonds()
        bondMask[5] = False
        bondMask[6] = False
        bondMask[7] = False

        visitor = helium.DFSDebugVisitor()
        helium.depth_first_search_mask(mol, mol.atom(0), visitor, atomMask, bondMask)
        self.compare_file('dfs9.log', visitor.output)

    def test_dfs_10(self):
        mol = helium.Molecule()
        SMILES.read('C(C)N', mol)

        visitor = helium.DFSDebugVisitor()
        helium.ordered_depth_first_search(mol, [0, 2, 1], visitor)
        self.compare_file('dfs10.log', visitor.output)

    def test_dfs_11(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        visitor = helium.DFSDebugVisitor()
        helium.exhaustive_depth_first_search(mol, mol.atom(1), visitor)
        self.compare_file('dfs11.log', visitor.output)

    def test_atom_order_visitor(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        visitor = helium.DFSAtomOrderVisitor()
        helium.depth_first_search(mol, visitor)
        self.assertListEqual([0, 1, 2], visitor.atoms)

    def test_bond_order_visitor(self):
        mol = helium.Molecule()
        SMILES.read('CCC', mol)

        visitor = helium.DFSBondOrderVisitor()
        helium.depth_first_search(mol, visitor)
        self.assertListEqual([0, 1], visitor.bonds)

    def test_closure_recorder_visitor(self):
        mol = helium.Molecule()
        SMILES.read('C1CCCCC1', mol)

        visitor = helium.DFSClosureRecorderVisitor()
        helium.depth_first_search(mol, visitor)
        self.assertListEqual([5], visitor.back_bonds)

    def test_custom_visitor(self):
        mol = helium.Molecule()
        SMILES.read('C1CCC(CC)CC1', mol)

        visitor = DFSTestVisitor()
        helium.depth_first_search(mol, visitor)
        self.compare_file('dfs2.log', visitor.output)


if __name__ == '__main__':
    unittest.main()
