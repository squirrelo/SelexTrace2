from unittest import TestCase, main
from os import tmpfile, remove
from multiprocessing import Manager

from selextrace.ctilib import (fold_clusters, SeqStructure, build_reference,
                               group_to_reference, group, group_by_seqstruct,
                               final_fold, create_group_output, make_r2r,
                               score_local_rnaforester, score_multi_forester,
                               group_by_forester, run_infernal)


class MainTests(TestCase):
    def setUp(self):
        self.del_files = []
        self.fastaseqs = [("seq1 count:11",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACAC"
                           "AUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGGGCAU"),
                          ("seq10 count:2",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAGCGCGUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGGGCAU"),
                          ("seq11 count:2",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAA"),
                          ("seq12 count:2",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACCCGGGGGCAU"),
                          ("seq2 count:26",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGCGCGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU"),
                          ("seq3 count:90",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU"),
                          ("seq4 count:17",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGCGCAU"),
                          ("seq5 count:7",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGCCCUAACACGGGGGCAU"),
                          ("seq6 count:4",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGCCCGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU"),
                          ("seq7 count:3",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGCCCUAACACGGGCGCAU"),
                          ("seq8 count:2",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGCACAU"),
                          ("seq9 count:2",
                           "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGAGUGGAUAUGGUACGC"
                           "AUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU")]
        self.fastastruct = ".....((((((..((((((((...((....(((.....)))..))...))))))....))..)).))))....((((((......))))))."
    def tearDown(self):
        for f in self.del_files:
            remove(f)

    def test_fold_clusters(self):
        """ Makes sure fold_clusters function returns proper folded structure
        """
        cfo = open("cluster_structs.fasta", 'w')
        cfo.close()
        manager = Manager()
        lock = manager.Lock()
        fold_clusters(lock, "cluster_1", self.fastaseqs, "./")
        with open("cluster_structs.fasta") as obsin:
            obs_fold = obsin.read()
        self.del_files.append("cluster_structs.fasta")
        expected = ">cluster_1\n%s\n" % self.fastastruct
        for header, seq in self.fastaseqs:
            expected = ''.join([expected, ">%s\n%s\n" % (header, seq)])
        self.assertEqual(obs_fold, expected)

    def test_build_reference_lengths(self):
        """Tests build_reference creates proper length ref and nonref lists"""
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        obs1, obs2 = build_reference(items, 3)
        self.assertEqual(len(obs1), 3)
        self.assertEqual(len(obs2), 7)

    def test_build_reference_items(self):
        """Tests whether items in build_reference ref and nonref lists are
        correct """
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        obs1, obs2 = build_reference(items, 3)
        self.assertEqual(len(obs1), 3)
        self.assertEqual(len(obs2), 7)
        #check that the ref and nonref are same
        for num in obs1:
            if num not in items or num in obs2:
                raise AssertionError(str(num) + " not expected in ref!")
        for num in obs2:
            if num not in items or num in obs1:
                raise AssertionError(str(num) + " not expected in nonref!")
        for num in items:
            if num not in obs1 and num not in obs2:
                raise AssertionError(str(num) + " not observed!")

    def test_build_reference_dupes(self):
        """Tests whether duplicates introduced by build_reference"""
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        obs1, obs2 = build_reference(items, 3)
        self.assertEqual(len(obs1), 3)
        self.assertEqual(len(obs2), 7)
        #check that the ref and nonref are same
        finals = set([])
        for item in obs1:
            if item in finals:
                raise AssertionError("Duplicate in reference!")
            finals.add(item)
        for item in obs2:
            if item in finals:
                raise AssertionError("Duplicate in nonreference!")
            finals.add(item)

    def test_group_to_reference(self):
        """Tests groups created by group_to_reference"""
        pass


if __name__ == "__main__":
    main()
