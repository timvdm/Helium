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

        f = helium.BinaryInputFile()
        self.assertFalse(f)

        f = helium.BinaryInputFile()
        self.assertTrue(f.open(self.datadir + '/1K.hel'))
        self.assertTrue(f)

        f = helium.BinaryInputFile()
        self.assertRaises(RuntimeError, f.open, 'foo')
        self.assertFalse(f)

        f = helium.BinaryInputFile(self.datadir + '/1K.hel')
        self.assertTrue(f)

        header = f.header()
        self.assertTrue(isinstance(header, str))

        f.close()
        self.assertFalse(f)


if __name__ == '__main__':
    unittest.main()
