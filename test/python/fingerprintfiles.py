import helium
import unittest
import os

class TestBinaryInputFile(unittest.TestCase):

    datadir = os.getenv('HEDATADIR', os.getenv('HOME') + '/Helium/data')

    def test_data_dir(self):
        self.assertTrue(os.path.isdir(self.datadir))

    def test_in_memory_row_major_fingerprint_storage(self):
        storage = helium.InMemoryRowMajorFingerprintStorage()
        self.assertFalse(storage.load('foo'))

        storage.load(self.datadir + '/1K_fp_path_7_1024_row.hel')
        self.assertTrue(isinstance(storage.header(), str))
        self.assertEqual(1024, storage.numBits())
        self.assertEqual(1000, storage.numFingerprints())
        self.assertTrue(isinstance(storage.fingerprint(0), helium.Fingerprint))

    def test_in_memory_column_major_fingerprint_storage(self):
        storage = helium.InMemoryColumnMajorFingerprintStorage()
        self.assertFalse(storage.load('foo'))

        storage.load(self.datadir + '/1K_fp_path_7_1024_col.hel')
        self.assertTrue(isinstance(storage.header(), str))
        self.assertEqual(1024, storage.numBits())
        self.assertEqual(1000, storage.numFingerprints())
        self.assertTrue(isinstance(storage.bit(0), helium.Fingerprint))

    def test_row_major_fingerprint_output_file(self):
        f = helium.RowMajorFingerprintOutputFile('tmp', 1021)
        fp = helium.Fingerprint(16)
        f.writeFingerprint(fp)
        f.writeFingerprint(fp)
        header = '{ "filetype": "fingerprints", "num_bits": 1021, "num_fingerprints": 2, "order": "row-major" }'
        f.writeHeader(header)

        storage = helium.InMemoryRowMajorFingerprintStorage()
        storage.load('tmp')
        self.assertEqual(header, storage.header())
        self.assertEqual(1021, storage.numBits())
        self.assertEqual(2, storage.numFingerprints())

    def test_column_major_fingerprint_output_file(self):
        f = helium.ColumnMajorFingerprintOutputFile('tmp', 1021, 2)
        fp = helium.Fingerprint(16)
        f.writeFingerprint(fp)
        f.writeFingerprint(fp)
        header = '{ "filetype": "fingerprints", "num_bits": 1021, "num_fingerprints": 2, "order": "column-major" }'
        f.writeHeader(header)

        storage = helium.InMemoryColumnMajorFingerprintStorage()
        storage.load('tmp')
        self.assertEqual(header, storage.header())
        self.assertEqual(1021, storage.numBits())
        self.assertEqual(2, storage.numFingerprints())



if __name__ == '__main__':
    unittest.main()
