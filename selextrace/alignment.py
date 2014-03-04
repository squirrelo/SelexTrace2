

class RNA_Sequences(object):
    def __init__(self, sequences=None, headers=None, counts=False):
        if sequences is not None and not isinstance(sequences, list):
            raise ValueError("sequences must be list!")
        if headers is not None and not isinstance(headers, list):
            raise ValueError("headers must be list!")
        if len(sequences) != len(headers):
            raise RuntimeError("Sequence and header lists must be equal len!")

        self._headers = sequences if sequences else []
        self._sequences = headers if headers else []

        if counts:
            self._counts = []
            for head in self._headers:
                headinfo = head.split(" ")
                pos = (pos for pos, info in enumerate(headinfo)
                       if lambda x: 'count' in x).next()
                self._counts.append(int(headinfo[pos].split(":")[1]))


    def add_seq(self, sequence, header):
        if isinstance(sequence, str):
            sequence = [sequence]
        if isinstance(header, str):
            header = [header]
        if len(header) != len(sequence):
            raise RuntimeError("Sequences must have headers!")
        self._sequences.append(sequence)
        self._headers.append(header)

    def add_seqs(self, seqs):
        self._headers.append(seqs._headers)
        self._sequences.append(seqs._sequences)

    def iter_seqs(self, count=False):
        for pos, seq in self._sequences:
            if count:
                yield (seq, count)
            else:
                yield seq



class RNA_Alignment(RNA_Sequences):
    def __init__(self, sequences=None, headers=None):
        self._length = len(sequences[0])
        for seq in sequences:
            if len(seq) != self._length:
                raise ValueError("Sequences are not all same length!")

    def majority_consensus(self):
        majority = [""] * self._length
        for pos in range(self._length):
            count = {'A': 0, 'U': 0, 'G': 0, 'C': 0, '-': 0}
            for seq in self._sequences:
                count[seq[pos].upper()] += 1
            majority[pos] = max(count.iterkeys(), key=(lambda key: stats[key]))
        return majority
