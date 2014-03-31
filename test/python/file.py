import helium
import unittest
import os

class TestBinaryInputFile(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def test_file(self):
        f = helium.BinaryInputFile()

        self.assertFalse(f)




if __name__ == '__main__':
    unittest.main()
