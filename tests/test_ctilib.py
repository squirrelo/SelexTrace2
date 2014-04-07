from unittest import TestCase, main
from os import remove
from multiprocessing import Manager

from selextrace.ctilib import (fold_clusters, SeqStructure, build_reference,
                               group_to_reference, group, group_by_seqstruct,
                               make_r2r, run_infernal, create_seqstructs)


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
        self.seqstruct = [SeqStructure("(((....)))...((((.((.(((((.........((((((((...(((...((......))...))).))))))))..))))).))))))...(((((((....)))))))..", "GACUUCGGUCCAAGCUAAUGCACUCACACAGACUCGUGGAUAUGGCACGCUACCUACGCUGGGACCGUAAUGUCCAUUAUGGGUUCAUAGCAACAUCGGGCUUCGGUCCGGUUC", "cluster_338"),
                          SeqStructure("....((((.(((((((.((((((((........((.((.......)).)).......(((((((((.......))))).)))).))))))))..))))).)).))))...", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGCACGCAUUGAAUUGUUGGGCACCGCAAAUGCCCUAACACGGGUGCAUCGGGCUUCGGUCCGGUUC", "cluster_340"),
                          SeqStructure(".....(((.(((((((..(((((.(((......(((((.......))))).......(((((((((.......))))).))))))))))))...))))).)).)))....", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUACGGCACGCAUUGAAUUGUUGGACACCUUAAAUGUCCUAACACGGGUGCAACGGGCUUCGGUCCGGUUC", "cluster_343"),
                          SeqStructure("....((((.(((((((.((((((.(((...(((.....)))................(((((((((.......))))).)))))))))))))..))))).)).))))...", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGACACCUUAAAUGUCCUAACACGGGUGCAUCGGGCUUCGGUCCGGUUC", "cluster_344"),
                          SeqStructure(".....(((.(((((((.((((((((((......(((((.......)))))...........(((((.......)))))...)).))))))))..))))).)).)))....", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGCUGGACACCUUAAAUGUCCUAACACGGGUGCAUCGGGCUUCGGUCCGGUUC", "cluster_345"),
                          SeqStructure("....((((.(((((((.((((((.(((...(((.....)))(((.....))).....(((((((.((.......))))))))))))))))))..))))).)).))))...", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACACAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGUGCAUCGGGCUUCGGUCCGGUUC", "cluster_339"),
                          SeqStructure("(((....)))...(((.(((.((((...(((..(((((......)))))..(((.....))).............)))....)))).))).)))..(((((((....)))))))..", "GACUUCGGUCCAAGCUAAUGCACUCUUAACGCCGCGUGGAUAUGCACGCAACCGUGAAUCGGACACCGUAAAUUCCGUAAGUGGGUACAUCAGCAAAUCGGGCUUCGGUCCGGUUC", "cluster_342"),
                          SeqStructure("....((((.(((((((.((((.(((.......(((((.......)))))........(((((((((.......))))).)))).))).))))..))))).)).))))...", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGUGUGGAUAUGGCACGCAUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGUGCAUCGGGCUUCGGUCCGGUUC", "cluster_337"),
                          SeqStructure("(((....)))...((((((.((((((........(((((((...((((...((......))...)))).)))))))..))))))))))))...(((((((....)))))))..", "GACUUCGGUCCAAGCUAAUGCACUCCCAUUUUCCGUGGAUAUGGUACGCUACCUACGUUGGGACCGUAAUGUCCACUAGGGGUGAUUAGCAAAAUCGGGCUUCGGUCCGGUUC", "cluster_346"),
                          SeqStructure("(((....)))...((.(((((((((..........((((((((...(((...((......))...))).))))))))...)))))))))))(((((((....)))))))..", "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGCACGCUACCUAAUACGGGACCGUAAUGUCCAUUAUGGGUGCAUUGCAUCGGGCUUCGGUCCGGUUC", "cluster_347")]

    def tearDown(self):
        for f in self.del_files:
            remove(f)

    def test_fold_clusters(self):
        """ Makes sure fold_clusters function returns proper folded structure
        """
        cfo = open("./cstest.fasta", 'w')
        cfo.close()
        manager = Manager()
        lock = manager.Lock()
        fold_clusters(lock, "cluster_1", self.fastaseqs, "./cstest.fasta")
        with open("./cstest.fasta") as obsin:
            obs_fold = obsin.read()
        self.del_files.append("./cstest.fasta")
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

    def test_create_seqstructs(self):
        """tests that create_seqstructs makes SeqStructure objects list"""
        with open("./support_files/cs.fasta") as fin:
            obs = create_seqstructs(fin, 10)
            self.assertEqual(obs, self.seqstruct)

    def test_create_seqstructs_wrong_numclusts(self):
        """tests that create_seqstructs thows error if incorrect number of
        SeqStructure objects"""
        with open("./support_files/cs.fasta") as fin:
            self.assertRaises(AssertionError, create_seqstructs, fin, 300)

    def test_group(self):
        """Tests groups created by group"""
        obs_group, obs_nogroup = group(self.seqstruct, 0.75)
        exp_group = {'cluster_337': ['cluster_343', 'cluster_345',
                                     'cluster_339'],
                     'cluster_347': ['cluster_338'],
                     'cluster_344': ['cluster_340']}
        exp_nogroup = [self.seqstruct[6], self.seqstruct[8]]

        self.assertEqual(obs_group, exp_group)
        self.assertEqual(obs_nogroup, exp_nogroup)

    def test_group_with_ref(self):
        obs_group, obs_nogroup = group(self.seqstruct[3:], 0.75,
                                       self.seqstruct[:3])
        exp_group = {'cluster_343': ['cluster_345', 'cluster_337'],
                     'cluster_338': [],
                     'cluster_340': ['cluster_344', 'cluster_339']}
        exp_nogroup = [self.seqstruct[6], self.seqstruct[8], self.seqstruct[9]]

        self.assertEqual(obs_group, exp_group)
        self.assertEqual(obs_nogroup, exp_nogroup)

    def test_group_to_reference(self):
        obs_group, obs_nogroup = group_to_reference(self.seqstruct[:3],
                                                    self.seqstruct[3:],
                                                    0.75, cpus=3)
        exp_group = {'cluster_343': ['cluster_337', 'cluster_345'],
                     'cluster_338': [],
                     'cluster_340': ['cluster_339', 'cluster_344']}
        exp_nogroup = [self.seqstruct[6], self.seqstruct[8], self.seqstruct[9]]
        #due to nature of multiprocess, list can be in any order
        #need to sorth then check to make sure correct
        for key in obs_group:
            obs_group[key].sort()
        obs_nogroup.sort(key=lambda x: x.name)
        self.assertEqual(obs_nogroup, exp_nogroup)
        self.assertEqual(obs_group, exp_group)

if __name__ == "__main__":
    main()
