from cogent import LoadSeqs
from cogent.app.uclust import Uclust, clusters_from_uc_file
from .qiime import local_align_primer_seq


def write_fasta_list(lst, filename):
    '''writes formatted list [(header,sequence)] to fasta file filename'''
    fileout = open(filename, 'w')
    for seq in lst:
        fileout.write('>%s\n%s\n' % (seq[0], seq[1]))
    fileout.close()


def write_fasta_dict(dct, filename):
    '''writes formatted dict {header: sequence} to fasta file filename'''
    fileout = open(filename, 'w')
    for header in dct:
        fileout.write('>%s\n%s\n' % (header, dct[header]))
    fileout.close()


def strip_primer(seqs, primer, maxmismatch=0, keep_primer=False):
    '''strips 3 prime primer from sequences in fasta file and returns arrays
        for stripped and not stripped sequences'''
    nostripped = []
    stripped = []
    pri = primer.upper()
    for head, seq in seqs:
        RNA = False
        seq = seq.upper()
        if 'U' in seq:
            seq = seq.replace('U', 'T')
            RNA = True
        #code adapted from truncate_reverse_primers.py in qiime
        rev_primer_mm, rev_primer_index = local_align_primer_seq(pri, seq)
        if rev_primer_mm > maxmismatch:
            nostripped.append((head, seq))
            continue
        if keep_primer:
            seqnew = seq[:rev_primer_index + len(primer)]
        else:
            seqnew = seq[:rev_primer_index]
        if RNA:
            seqnew = seqnew.replace('T', 'U')
        stripped.append((head, seqnew))
    #end for
    return stripped, nostripped


def rem_N_short(seqs, minlen=1):
    '''Takes in a [(header, seq)] formatted list of sequences and returns list
    with sequences containing Ns or shorter than minlen removed'''
    rem = []
    for i, seq in enumerate(seqs):
        if "N" in seq[1].upper() or len(seq[1]) < minlen:
            rem.append(i)
    rem.sort(reverse=True)
    for i in rem:
        seqs.pop(i)
    return seqs


def cluster_seqs(seqspath, simm, folderout='/tmp', gapopen=None, gapext=None):
    if folderout[-1] != "/":
        folderout += "/"

    params = {
        '--usersort': False,
        '--id': float(simm),
        '--maxaccepts': 20,
        '--maxrejects': 500,
        '--stepwords': 20,
        '--hsp': 0
    }
    if gapopen is not None:
        params['--gapopen'] = gapopen
    if gapext is not None:
        params['--gapext'] = gapext
    uclust = Uclust(params, WorkingDir='/tmp')
    input_data = {
        '--input': seqspath,
        '--uc': folderout + "clusters.uc",
        '--log': folderout + "clusters.log"
    }
    result = uclust(input_data)
    clusters, failures, new_seeds = clusters_from_uc_file(result['ClusterFile'])

    seqs = LoadSeqs(seqspath, aligned=False)
    convheader = {}
    clusterseqs = {}
    #create dictinary to convert shortened headers to full headers
    for header in seqs.getSeqNames():
        convheader[header.split()[0]] = header
    #match headers in each cluster to seqs to create cluster tuples list
    for num, cluster in enumerate(clusters):
        clusterseqs["cluster_" + str(num)] = []
        for header in clusters[cluster]:
            clusterseqs["cluster_" + str(num)].append((convheader[header],
                                        seqs.getSeq(convheader[header])))

    return clusterseqs

def count_seqs(headers, field="count", divider=":"):
    """Returns the total sequence count of deduplicated sequences"""
    count = 0
    if not isinstance(headers, list):
            headers = [headers]
    for header in headers:
        head = header.split()
        hold = next((val for val in head if field in val), None)
        if hold is None:
            raise RuntimeError("%s not found in header %s" % (field, header))
        count += int(hold.split(divider)[1])
    return count

def remove_duplicates(seqsin):
    '''Takes in LoadSeqs loadable sequences, removes duplicate sequences
    and returns a list of unique sequence tuples, formated (sequence, count)
    sorted most abundant to least abundant'''

    parsable_seqs = LoadSeqs(data=seqsin, aligned=False)
    uniques = {}
    for header, seq in parsable_seqs.items():
        seq = str(seq)
        if seq in uniques:
            uniques[seq] += 1
        else:
            uniques[seq] = 1
    uniques = uniques.items()
    uniques.sort(key=lambda x: x[1], reverse=True)
    return uniques
