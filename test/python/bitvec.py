import helium
import unittest

class TestBitvec(unittest.TestCase):

    def test_Fingerprint(self):
        fp = helium.Fingerprint(8)
        self.assertEqual(8, fp.numWords)

        fp = helium.Fingerprint(1)
        self.assertEqual('0000000000000000', fp.hex())
        self.assertEqual('0000000000000000 0000000000000000 0000000000000000 0000000000000000', fp.bin())
        self.assertEqual('0000000000000000000000000000000000000000000000000000000000000000', fp.bin(False))

        fp = helium.bitvec_from_binary('0101010101010101010010101010')
        self.assertEqual('0101010101010101010010101010000000000000000000000000000000000000', fp.bin(False))
        self.assertEqual(1, fp.numWords)

        fp = helium.bitvec_from_hex('0123456789abcdefABCDEF')
        self.assertEqual('0123456789abcdefabcdef0000000000', fp.hex())
        self.assertEqual(2, fp.numWords)

        fp = helium.Fingerprint(1)
        fp[4] = 1
        self.assertEqual(1, fp[4])
        self.assertEqual(1, fp.count())
        fp[4] = 0
        self.assertEqual(0, fp.count())
        fp[4] = 1
        fp[8] = 1
        fp[12] = 1
        fp[40] = 1
        self.assertEqual(4, fp.count())
        self.assertEqual(2, fp.count(6, 30))
        fp.zero()
        self.assertEqual(0, fp.count())

    def test_bitvec(self):
        self.assertEqual(64, helium.BitsPerWord)
        self.assertEqual(1, helium.bitvec_num_words_for_bits(42))
        self.assertEqual(2, helium.bitvec_num_words_for_bits(70))

        fp = helium.Fingerprint(1)

        fp2 = helium.bitvec_copy(fp)
        helium.bitvec_set(4, fp2)
        self.assertTrue(helium.bitvec_get(4, fp2))
        self.assertFalse(helium.bitvec_get(4, fp))
        helium.bitvec_zero(fp2)
        self.assertFalse(helium.bitvec_get(4, fp2))

        helium.bitvec_set(4, fp)
        helium.bitvec_reset(4, fp)
        self.assertFalse(helium.bitvec_get(4, fp))

        helium.bitvec_set(4, fp)
        helium.bitvec_set(8, fp)
        helium.bitvec_set(4, fp2)
        helium.bitvec_set(8, fp2)
        helium.bitvec_set(12, fp2)
        self.assertTrue(helium.bitvec_is_subset_superset(fp, fp2))
        self.assertFalse(helium.bitvec_is_subset_superset(fp2, fp))
        self.assertEqual(2, helium.bitvec_count(fp))
        self.assertEqual(3, helium.bitvec_count(fp2))
        helium.bitvec_set(40, fp2)
        self.assertEqual(3, helium.bitvec_count(fp2, 0, 32))
        self.assertEqual(1, helium.bitvec_count(fp2, 32, 64))

        helium.bitvec_union_count(fp, fp2)
        helium.bitvec_tanimoto(fp, fp2)
        helium.bitvec_tanimoto(fp, fp2, 0, 1)
        helium.bitvec_cosine(fp, fp2)
        helium.bitvec_cosine(fp, fp2, 0, 1)
        helium.bitvec_hamming(fp, fp2)
        helium.bitvec_hamming(fp, fp2, 0, 1)
        helium.bitvec_russell_rao(fp, fp2)
        helium.bitvec_forbes(fp, fp2)
        helium.bitvec_forbes(fp, fp2, 0, 1)




if __name__ == '__main__':
    unittest.main()
