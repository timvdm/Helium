import helium
import unittest

class TestElement(unittest.TestCase):

    def test_symbol(self):
        self.assertEqual('C', helium.Element.symbol(6))
        self.assertEqual('O', helium.Element.symbol(8))
        self.assertEqual('S', helium.Element.symbol(16))

    def test_element(self):
        for i in range(112):
            symbol = helium.Element.symbol(i + 1)
            self.assertEqual(i + 1, helium.Element.element(symbol))

    def test_average_mass(self):
        self.assertEqual(12, helium.Element.averageMass(6))
        self.assertEqual(16, helium.Element.averageMass(8))
        self.assertEqual(32, helium.Element.averageMass(16))

    def test_addHydrogens(self):
        self.assertTrue(helium.Element.addHydrogens(6))
        self.assertTrue(helium.Element.addHydrogens(8))
        self.assertFalse(helium.Element.addHydrogens(21))

    def test_valence(self):
        self.assertEqual(4, helium.Element.valence(6, 0, 0))
        self.assertEqual(4, helium.Element.valence(6, 0, 3))
        self.assertEqual(3, helium.Element.valence(6, -1, 3))

if __name__ == '__main__':
    unittest.main()
