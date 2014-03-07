from unittest import TestCase, main
from selextrace.stutils import (write_fasta_list, write_fasta_dict,
                                strip_primer, rem_N_short, cluster_seqs,
                                count_seqs, remove_duplicates)


class MainTests(TestCase):
    def test_write_fasta_list(self):
        pass

    def test_write_fasta_dict(self):
        pass

    def test_strip_primer_keep(self):
        pass

    def test_strip_primer_nokeep(self):
        pass

    def test_rem_N_short(self):
        pass

    def test_rem_N_short_nolen(self):
        pass

    def test_cluster_seqs(self):
        pass

    def test_count_seqs(self):
        pass

    def test_remove_duplicates(self):
        testin = [
            ("test1", "AAAAGGGCCCTTTAGCTAAA"),
            ("test2", "AAAGGGCCCTTTAAA"),
            ("test3", "AAAGGGCCCTTTAAA")]
        obs = remove_duplicates(testin)

        expected = [('AAAGGGCCCTTTAAA', 2),
                    ('AAAAGGGCCCTTTAGCTAAA', 1)]
        self.assertEqual(obs, expected)


if __name__ == "__main__":
    main()
