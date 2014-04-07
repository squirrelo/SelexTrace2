#!/usr/bin/env python
#file split_libraries.py

"""Code taken from split_libraries.py of qiime 1.8.0. All credits below left
    intact.
"""

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = [
    "Rob Knight",
    "Micah Hamady",
    "Greg Caporaso",
    "Kyle Bittinger",
    "Jesse Stombaugh",
    "William Walters",
    "Jens Reeder",
    "Emily TerAvest"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "rob@spot.colorado.edu, william.a.walters@colorado.edu"
from cogent import LoadSeqs, DNA
from cogent.align.align import make_dna_scoring_dict, local_pairwise
from cogent.core.moltype import IUPAC_DNA_ambiguities


def pair_hmm_align_unaligned_seqs(seqs, moltype=DNA, params={}):
    """
        Checks parameters for pairwise alignment, returns alignment.

        Code from Greg Caporaso.
    """

    seqs = LoadSeqs(data=seqs, moltype=moltype, aligned=False)
    try:
        s1, s2 = seqs.values()
    except ValueError:
        raise ValueError(
            "Pairwise aligning of seqs requires exactly two seqs.")

    try:
        gap_open = params['gap_open']
    except KeyError:
        gap_open = 5
    try:
        gap_extend = params['gap_extend']
    except KeyError:
        gap_extend = 2
    try:
        score_matrix = params['score_matrix']
    except KeyError:
        score_matrix = make_dna_scoring_dict(
            match=1, transition=-1, transversion=-1)

    return local_pairwise(s1, s2, score_matrix, gap_open, gap_extend)


def MatchScorerAmbigs(match, mismatch, matches=None):
    """ Alternative scorer factory for sw_align which allows match to ambiguous chars

    It allows for matching to ambiguous characters which is useful for
     primer/sequence matching. Not sure what should happen with gaps, but they
     shouldn't be passed to this function anyway. Currently a gap will only match
     a gap.

    match and mismatch should both be numbers. Typically, match should be
    positive and mismatch should be negative.

    Resulting function has signature f(x,y) -> number.

    Code original from Greg Caporaso
    """

    matches = matches or \
        {'A': {'A': None}, 'G': {'G': None}, 'C': {'C': None},
         'T': {'T': None}, '-': {'-': None}}
    for ambig, chars in IUPAC_DNA_ambiguities.items():
        try:
            matches[ambig].update({}.fromkeys(chars))
        except KeyError:
            matches[ambig] = {}.fromkeys(chars)

        for char in chars:
            try:
                matches[char].update({ambig: None})
            except KeyError:
                matches[char] = {ambig: None}

    def scorer(x, y):
        # need a better way to disallow unknown characters (could
        # try/except for a KeyError on the next step, but that would only
        # test one of the characters)
        if x not in matches or y not in matches:
            raise ValueError("Unknown character: %s or %s" % (x, y))
        if y in matches[x]:
            return match
        else:
            return mismatch
    return scorer

equality_scorer_ambigs = MatchScorerAmbigs(1, -1)


def local_align_primer_seq(primer, sequence, sw_scorer=equality_scorer_ambigs):
    """Perform local alignment of primer and sequence

        primer: Input primer sequence
        sequence: target sequence to test primer against

        Returns the number of mismatches,
         and the start position in sequence of the hit.

        Modified from code written by Greg Caporaso.
    """

    query_primer = primer

    query_sequence = str(sequence)

    # Get alignment object from primer, target sequence
    alignment = pair_hmm_align_unaligned_seqs([query_primer, query_sequence])

    # Extract sequence of primer, target site, may have gaps in insertions
    # or deletions have occurred.
    primer_hit = str(alignment.Seqs[0])
    target_hit = str(alignment.Seqs[1])

    # Count insertions and deletions
    insertions = primer_hit.count('-')
    deletions = target_hit.count('-')

    mismatches = 0
    for i in range(len(target_hit)):
        # using the scoring function to check for
        # matches, but might want to just access the dict
        if sw_scorer(target_hit[i], primer_hit[i]) == -1 and \
                target_hit[i] != '-' and primer_hit[i] != '-':
            mismatches += 1
    try:
        hit_start = query_sequence.index(target_hit.replace('-', ''))
    except ValueError:
        raise ValueError(
            'substring not found, query string %s, target_hit %s' %
            (query_sequence, target_hit))

    # sum total mismatches
    mismatch_count = insertions + deletions + mismatches

    return mismatch_count, hit_start
