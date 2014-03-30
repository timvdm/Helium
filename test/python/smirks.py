import helium
import unittest

SMILES = helium.Smiles()

class TestSmarts(unittest.TestCase):

    def test_valid_smirks(self):
        smirks = helium.Smirks()
        self.assertTrue(smirks.init('[C:1]', '[N:1]'))
        self.assertFalse(smirks.error())
        self.assertEqual(0, len(str(smirks.error())))

        self.assertTrue(smirks.init('[C:1]>>[N:1]'))
        self.assertFalse(smirks.error())
        self.assertEqual(0, len(str(smirks.error())))

    def test_invalid_smirks(self):
        smirks = helium.Smirks()
        self.assertFalse(smirks.init('fas', '[C]'))
        self.assertTrue(smirks.error().__nonzero__())
        self.assertNotEqual(0, len(str(smirks.error())))

    def test_fix(self):
        smirks = helium.Smirks()
        smirks.setFixMass(True)
        smirks.setFixHydrogens(False)

    def test_apply(self):
        smirks = helium.Smirks()
        self.assertTrue(smirks.init('[C:1]>>[N:1]'))
        mol = helium.Molecule()

        SMILES.read('C', mol)
        self.assertTrue(smirks.apply(mol, helium.RingSet(mol)))
        self.assertEqual('N', SMILES.write(mol))

        SMILES.read('O', mol)
        self.assertFalse(smirks.apply(mol, helium.RingSet(mol)))
        self.assertEqual('O', SMILES.write(mol))


    def test_requires_cycles(self):
        smirks = helium.Smirks()

        smirks.init('[C:1]', '[N:1]')
        self.assertFalse(smirks.requiresCycles())

        smirks.init('[CR:1]', '[N:1]')
        self.assertTrue(smirks.requiresCycles())

    def test_requires_explicit_hydrogens(self):
        smirks = helium.Smirks()

        smirks.init('[C:1]', '[N:1]')
        self.assertFalse(smirks.requiresExplicitHydrogens())

        smirks.init('[C:1][H:2]', '[N:1][H:2]')
        self.assertTrue(smirks.requiresExplicitHydrogens())

if __name__ == '__main__':
    unittest.main()
