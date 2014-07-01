import helium
import unittest
import os

class TestBinaryInputFile(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME', 'c:/') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def test_fps(self):
        pass

if __name__ == '__main__':
    unittest.main()
