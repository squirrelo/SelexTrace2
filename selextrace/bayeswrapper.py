from cogent import LoadSeqs, RNA
from cogent.app.muscle_v38 import align_unaligned_seqs
from Bayes.bayes import BayesCalculation


class BayesInputWrapper:
    def __init__(self, labels, seqs, temp='37'):
        self.name = 'name'
        self.temperature = temp
        self.primer3 = ''
        self.primer5 = ''
        self.labels = labels
        self.sequences = seqs
        self.mappings = []
        self.id = '12345'


def bayesfold(seqsin, temperature=37, params=None, align=True):
    '''Takes in sequences in LoadSeqs readable format and returns
    most likely structure from bayesfold.'''
    try:
        if params is None:
            params = {}
        if not align:
            aln = LoadSeqs(data=seqsin, moltype=RNA, aligned=True)
        else:
            aln = align_unaligned_seqs(seqsin, RNA, params=params)
        bayesinput = BayesInputWrapper(aln.getSeqNames(),
                                       map(str, aln.iterSeqs()),
                                       str(temperature))
        bayescalc = BayesCalculation(bayesinput)
        bayescalc.run()
        struct = str(bayescalc.Alignment.Structures).split()[1]
        del bayescalc
        del bayesinput
        return aln, struct
    except Exception, e:
        print "BAYESFOLD ERROR: ", e
        raise RuntimeError("BAYESFOLD ERROR: ", e)
