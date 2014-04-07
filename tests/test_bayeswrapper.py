from unittest import TestCase, main

from cogent import LoadSeqs, RNA

from selextrace.bayeswrapper import bayesfold


class MainTests(TestCase):
    def setUp(self):
        self.seqs = {"seq1":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACACAUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGGGCAU",
            "seq2":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGCGCGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "seq3":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "seq4":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGCGCAU",
            "seq5":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGCCCUAACACGGGGGCAU",
            "seq6":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGCCCGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "seq7":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGCCCUAACACGGGCGCAU",
            "seq8":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGCACAU",
            "seq9":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGAGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "seq10":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAGCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGGGCAU",
            "seq11":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAA",
            "seq12":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACCCGGGGGCAU"}

        self.badseqs = [
            ("seq1", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGG"),
            ("seq2", "GACUUCGGUCCAAGCUAAC")]

        self.aln = LoadSeqs(data=self.seqs, moltype=RNA)

        self.struct = ".....((((((..((((((((...((....(((.....)))..))...))))))....))..)).))))....((((((......))))))."

    def test_bayesfold(self):
        """test bayesfold with default paramters"""
        obs_aln, obs_struct = bayesfold(self.seqs)
        self.assertEqual(self.struct, obs_struct)
        self.assertEqual(self.aln, obs_aln)

    def test_bayesfold_params(self):
        """test bayesfold with passed alignment parameters"""
        params = {"-diags": True, "-maxiters": 5}
        obs_aln, obs_struct = bayesfold(self.seqs, params=params)
        self.assertEqual(self.struct, obs_struct)
        self.assertEqual(self.aln, obs_aln)

    def test_bayeswrapper_noalign(self):
        """Test bayeswrapper with previously aligned sequences"""
        obs_aln, obs_struct = bayesfold(self.aln, align=False)
        self.assertEqual(self.struct, obs_struct)
        self.assertEqual(self.aln, obs_aln)

    def test_bayeswrapper_noalign_badseqs(self):
        """Test bayeswrapper throws error when unaligned sequences passed as
        aligned"""
        self.assertRaises(RuntimeError, bayesfold, self.badseqs, align=False)


if __name__ == "__main__":
    main()
