from os.path import exists, join
from os import mkdir, walk
from subprocess import Popen, PIPE
from math import ceil
from random import shuffle
from multiprocessing import Pool
from traceback import format_exc
from collections import defaultdict

from cogent import LoadSeqs, RNA
from cogent.app.infernal_v11 import (cmsearch_from_file, calibrate_file)
from cogent.format.stockholm import stockholm_from_alignment
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.muscle_v38 import align_unaligned_seqs
from skbio.core.sequence import RNASequence
from skbio.core.alignment import StockholmAlignment
from numpy import empty, savetxt, ones, identity, int16
import align

from bayeswrapper import bayesfold
from selextrace.stutils import count_seqs


def fold_clusters(lock, cluster, seqs, otufile):
    '''Function for multithreading.
    Computes structure for a cluster and writes it to file'''
    aln, struct = bayesfold(seqs, params={"-diags": True})
    # write structure out to file
    try:
        lock.acquire()
        cfo = open(otufile, 'a')
        cfo.write(">%s\n%s\n" % (cluster, struct))
        cfo.write(aln.toFasta() + "\n")
        cfo.close()
        lock.release()
    except Exception:
        lock.release()


def write_clusters(clusters, cout):
    """ writes out cluster fasta file
    INPUT
    -----
    clusters: dict of list of tuples
        holds all clusters found
    cout: open file object
        open file to write into
    """
    hold = clusters.keys()
    hold.sort()
    cout.write(str(len(clusters)) + "\n")
    for cluster in hold:
        cout.write(">%s\n%s\n" % (cluster, cluster))
        for header, seq in clusters[cluster]:
            cout.write(">%s\n%s\n" % (header, seq))


def read_clusters(cfo):
    numclusts = int(cfo.readline().strip())
    currclust = ""
    clusters = defaultdict(list)
    for header, seq in MinimalFastaParser(cfo):
        if "cluster" in header:
            currclust = header
        else:
            clusters[currclust].append((header, seq))

    return clusters, numclusts


def create_seqstructs(cfo, numclusts):
    seqstructs = []
    # read in first cluster and struct
    currclust = cfo.readline().strip(">").strip()
    struct = cfo.readline().strip()
    seqs = []
    for header, seq in MinimalFastaParser(cfo):
        if "cluster_" in header:
            aln = LoadSeqs(data=seqs, moltype=RNA)
            seqstructs.append(SeqStructure(struct,
                                           ''.join(aln.majorityConsensus()),
                                           currclust))
            # move on to next structgroup
            struct = seq
            seqs = []
            currclust = header
        else:
            seqs.append((header, seq))
    aln = LoadSeqs(data=seqs, moltype=RNA)
    seqstructs.append(SeqStructure(struct, ''.join(aln.majorityConsensus()),
                                   currclust))
    if len(seqstructs) != numclusts:
        raise AssertionError("%i structures, %i clusters. Not all clusters "
                             "folded!" % (len(seqstructs), numclusts))
    return seqstructs


class SeqStructure(object):
    __slots__ = ("pairs", "struct", "seq", "name", "_seqmap")

    def __init__(self, structure, seq, name):
        if len(seq) != len(structure):
            raise ValueError("sequence must be same length as structure!")

        self.struct = structure
        self.seq = seq
        self.name = name
        self._seqmap = self._create_seqmap()

        self.pairs = ('AU', 'UA', 'GC', 'CG', 'GU', 'UG')

    def __len__(self):
        return len(self.struct)

    def __eq__(self, other):
        if self.struct != other.struct:
            return False
        if self.seq != other.seq:
            return False
        if self.name != other.name:
            return False
        return True

    def __ne__(self, other):
        if self.struct == other.struct:
            return False
        if self.seq == other.seq:
            return False
        if self.name == other.name:
            return False
        return True

    def _create_seqmap(self):
        seqmap = []
        stack = []
        for pos, symbol in enumerate(self.struct):
            if symbol == "(":
                stack.append(pos)
            elif symbol == ")":
                p1 = stack.pop()
                pair = ''.join((self.seq[p1], self.seq[pos]))
                seqmap.append((p1, pos, pair))
        return tuple(seqmap)

    def score_seq(self, seq):
        """Scores sequence based on how well it fits into given seq/struct map.
            INPUT:
            seq: RNA sequence to thread into structure
            OUTPUT
            Score, where pair is given 2 points if matches pair in structure
            and one point if it can base pair at all. This divided by seqmap
            length to get final score between 0 and 2.
        """
        maxscore = 0
        score = 0
        start = 0
        end = len(seq)
        if self._seqmap is None:
            raise RuntimeError("No seqmap exists!")
        if len(seq) == len(self.struct):
            # same length so trivial
            maxscore = self._eval_struct_seq(seq)
        elif len(seq) > len(self.struct):
            # sequence longer than struct so keep slicing sequence
            # evaluate all possible sequence fittings into the structure
            # return highest scoring one
            while end != len(seq):
                score = self._eval_struct_seq(seq[start:end])
                if maxscore < score:
                    maxscore = score
                start += 1
                end += 1
        else:
            # struct longer than sequence so walk seq down struct and compare
            # evaluate all possible sequence fittings into the structure
            # return highest scoring one
            while end != len(self.struct):
                score = self._eval_struct_seq(seq, seqoffset=start)
                if maxscore < score:
                    maxscore = score
                start += 1
                end += 1
        # return normalized score
        return float(maxscore) / len(self._seqmap)

    def _eval_struct_seq(self, seq, seqoffset=0):
        if self._seqmap is None:
            raise RuntimeError("There is no seqmap to evaluate!")
        score = 0
        for p1, p2, pair in self._seqmap:
            p2 = p2 - seqoffset
            # make sure we are in range of sequence given including offset
            if p1 < seqoffset or p2 >= len(seq):
                continue
            p1 = p1 - seqoffset
            # p1 & p2: positions of bases, pair: basepair at that position
            currpair = seq[p1] + seq[p2]
            currpair = currpair.upper()
            pair = pair.upper()
            if currpair in self.pairs:
                if currpair == pair:
                    # perfect pair match with bases scores 2 points
                    score += 2
                else:
                    # ability to basepair at all scores 1 point
                    score += 1
        return score


# object wrapper so create Popen object once: saves DAYS of overhead
class ScoreStructures(object):
    def __init__(self):
        self.p = Popen(["RNAdistance"], stdin=PIPE, stdout=PIPE)

    def __call__(self, struct1, struct2):
        strlen = len(struct1)
        if len(struct2) > strlen:
            strlen = len(struct2)
        self.p.stdin.write(''.join([struct1, "\n", struct2, "\n"]))
        self.p.stdin.flush()
        # lower number better for rnadistance program, so need to reverse for
        # our scoring algorithm
        score = strlen - float(self.p.stdout.readline().strip().split(":")[1])
        return score / strlen

    def end(self):
        self.p.stdin.write('@\n')
        self.p.stdin.flush()
        # force kill to make sure it closes
        try:
            self.p.kill()
        except:
            pass


def build_reference(keys, refsize):
    '''Creates a random list of references refsize long, returning rest as
       nonref
    '''
    shuffle(keys)
    # return reference, nonreference by slicing list
    return keys[:refsize], keys[refsize:]


def group_to_reference(reference, nonref, minscore, cpus=1):
    pool = Pool(processes=cpus)
    hold = []

    chunksize = int(ceil(len(nonref)/float(cpus)))
    if chunksize == 0:
        chunksize = 1
    for startpos in range(0, len(nonref), chunksize):
        # divide up nonref into chunks and align each chunk to reference seqs
        # final # chunks == number of cpus available
        endpos = startpos+chunksize
        pool.apply_async(func=group, args=(nonref[startpos:endpos], minscore,
                                           reference), callback=hold.append)
    pool.close()
    pool.join()
    # join all results to final dictionary and list to return
    groups = {r.name: [] for r in reference}
    nogroup = []
    for grouped, notgrouped in hold:
        nogroup.extend(notgrouped)
        for ref, g in grouped.iteritems():
            groups[ref].extend(g)
    return groups, nogroup


def group(nonref, minscore, ref=None):
    # takes in list of seqstructure objects for nonref and ref
    # if ref list is pased, know we are reference grouping

    denovo = False if ref else True
    nogroup = []
    grouped = defaultdict(list)

    # loop through all nonreference items
    for pos, currnonref in enumerate(nonref):
        if currnonref.name in grouped:
            continue
        seq1 = currnonref.seq.replace("-", "")
        bestref = None
        bestscore = minscore
        if denovo:
            ref = nonref[pos+1:]
        # compare to each reference item
        for refstruct in ref:
            seq2 = refstruct.seq.replace("-", "")
            # get alignment score and add to seq/struct score
            # aln, alnsc = nw_align(seq1, seq2, return_score=True)
            # score is normalized by dividing each score by sequence length
            # then adding. This should keep scores between zero and one
            align_score = nwalign_wrapper(seq1, seq2)
            struct_score = currnonref.score_seqstruct(refstruct.struct) / 2
            if (align_score >= minscore and struct_score >= minscore):
                score = (align_score + struct_score)
            else:
                # disregard as possible match
                score = -1
            if score > bestscore:
                bestscore = score
                bestref = refstruct.name

        if bestref:
            grouped[bestref].append(currnonref.name)
        elif currnonref.name not in grouped:
            nogroup.append(currnonref)

    # make sure all ref in the final grouping dictionary if applicable
    if ref is not None:
        for r in ref:
            if r.name not in grouped:
                grouped[r.name] = []
    return (grouped, nogroup)


NUCMATRIX = -ones((256, 256)) + 2 * identity(256)
NUCMATRIX = NUCMATRIX.astype(int16)


def nwalign_wrapper(seq1, seq2, matrix=NUCMATRIX):
    s1 = align.string_to_alignment(seq1)
    s2 = align.string_to_alignment(seq2)
    (score, a1, a2) = align.align(s1, s2, -1, -1, matrix)
    return float(score) / len(a1)


def group_by_seqstruct(grouping, structscore, setpercent=0.01, cpus=1):
        '''Does grouping by way of de-novo reference creation and clustering
            grouping - list of SeqStructure objects to be grouped
            structscore - maximum score to consider grouping structures
            specstructs - a list of a subset of structures in structgroups
                          to cluster (optional)
            setpercent - Allows manual percentage setting for ref structures
                           (default 1%  of dict)
        '''
        # fail if nothing to compare
        if len(grouping) < 1:
            raise ValueError("Must have at least one item to group!")
        # return the list directly if only one item (useful for breakout work)
        if len(grouping) == 1:
            return {grouping.name: []}
        # for speed, get 1% as initial clustering or user defined.
        # Need at least 5 structs though.
        finishlen = int(len(grouping) * setpercent)
        if finishlen == 0:
            finishlen = 1
        grouped = {}
        ungrouped = grouping
        #  keep refining while not at limit and are still grouping structs
        while len(ungrouped) > finishlen:
            ref, nonref = build_reference(ungrouped, finishlen)
            g, ungrouped = group_to_reference(ref, nonref, structscore, cpus)
            grouped.update(g)
        # do the last grouping
        g, ungrouped = group(ungrouped, structscore)
        grouped.update(g)
        # add ungroupable bit to end
        for ug in ungrouped:
            grouped[ug.name] = []
        return grouped


def align_order_seqs(seqs, params, outfolder, num, prefix="group_"):
    if exists("%s%s%i.fna" % (outfolder, prefix, num)):
        return
    try:
        aln = align_unaligned_seqs(seqs, RNA, params=params)
        aln.Names.sort(reverse=True, key=lambda c: count_seqs(c))
        with open("%s%s%i.fna" % (outfolder, prefix, num), 'w') as fout:
            fout.write(aln.toFasta() + "\n")
    except Exception as e:
        print("align_order_seqs ERROR: ", format_exc(e))


def create_final_output(groupfasta, basefolder, minseqs=1, cpus=1):
    '''Function for multithreading. Creates the final BayesFold alignment and
    writes to files, then r2r struct and infernal CM file'''
    # skip if already run and program just crashed or whatever
    currgroup = groupfasta.split("/")[-1].split(".")[0]
    currotufolder = basefolder + currgroup
    if exists(currotufolder):
        return

    # load seqs and make sure we have enough
    aln = LoadSeqs(groupfasta, moltype=RNA, aligned=True)
    count = count_seqs(aln.Names)
    if count < minseqs:
        return
    # get weights for each sequence. weight==count
    weights = []
    maxweight = 0
    for header in aln.Names:
        weight = count_seqs(header)
        if weight > maxweight:
            maxweight = weight
        weights.append(header.split()[0])
        weights.append(str(weight))

    # fold alignment with bayesfold
    aln, struct = bayesfold(aln, align=False)

    # write log information
    mkdir(currotufolder)
    with open(currotufolder + "/log.txt", 'w') as logout:
        logout.write(' '.join([currgroup, ":\n", str(count),
                     "sequences\n", str(aln.getNumSeqs()),
                              "unique sequences\nStructure: ", struct, "\n"]))
    # write out alignment and structure in fasta format
    with open(currotufolder + "/bayesfold-aln.fasta", 'w') as alnout:
        alnout.write(">SS_cons\n%s\n%s" % (struct, aln.toFasta()))

    # shave off info in header for stockholm
    aln = LoadSeqs(data=aln, moltype=RNA,
                   label_to_name=lambda x: x.split()[0])
    # create stockholm formatted alignment
    sto = stockholm_from_alignment(aln, GC_annotation={'SS_cons': struct})
    del aln
    # create standard weights for infernal
    infweights = ""
    for pos in range(0, len(weights), 2):
        infweights = ''.join([infweights, '# =GS %s WT %s\n' %
                             (weights[pos],
                              str(float(weights[pos+1]) / maxweight))])
    # create weights for r2r
    r2r_weights = "# =GF USE_THIS_WEIGHT_MAP " + ' '.join(weights)
    # create sto file with r2r and std weights
    sto = sto.split("\n")
    sto[-1] = infweights.strip()
    sto.append(r2r_weights)
    sto.append("//\n")
    stofile = currotufolder + "/bayesfold-aln.sto"
    with open(stofile, 'w') as alnout:
        alnout.write('\n'.join(sto))

    # make R2R secondary structure for alignment
    make_r2r(stofile, currotufolder, currgroup)
    # create CM file for infernal from group
    cmbuild_from_file(stofile, currotufolder + "/cmfile.cm",
                      params={'--wgiven': True})
    calibrate_cmfile(currotufolder + "/cmfile.cm", cpus=cpus)


def cmbuild_from_file(aln, cmout, params=None):
    if params is None:
        params = {}
    command = ["cmbuild"]
    for param in params:
        command.append(param)
        if not isinstance(params[param], bool):
            command.append(params[param])
    command.append(cmout)
    command.append(aln)
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    retcode = p.wait()
    if retcode != 0:
        raise RuntimeError("CM file build failed! %s" % p.stderr)


def calibrate_cmfile(cmfile, cpus=1):
    if not exists(cmfile):
        raise IOError("cmfile does not exist: %s" % cmfile)
    command = ["cmcalibrate", "--cpu", str(cpus), cmfile]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    retcode = p.wait()
    if retcode != 0:
        raise RuntimeError("CM file calibrate failed! %s" % ''.join(p.stderr))


def cmalign(inseqs, cmfile, outfile, cpus=1):
    '''runs cmalign for given information'''
    if not exists(cmfile):
        raise IOError("cmfile does not exist: %s" % cmfile)
    if not exists(inseqs):
        raise IOError("Sequence library does not exist: %s" % inseqs)
    command = ["cmalign", "--cpu", str(cpus), "-o", outfile, cmfile, inseqs]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    retcode = p.wait()
    if retcode != 0:
        raise RuntimeError("CM file calibrate failed! %s" % ''.join(p.stderr))


def make_r2r(insto, outfolder, group):
    '''generates R2R secondary structure pdf with default colorings'''
    try:
        outinfo = (outfolder, group)
        command = ["r2r", "--GSC-weighted-consensus", insto,
                   "%s/%s.sto" % outinfo, "3", "0.97", "0.9",
                   "0.75", "4", "0.97", "0.9", "0.75", "0.5", "0.1"]
        p = Popen(command)
        retcode = p.wait()
        if retcode != 0:
            raise RuntimeError("r2r: " + p.stderr)
        p = Popen(["r2r", "%s/%s.sto" % outinfo, "%s/%s.pdf" % outinfo],
                  stdout=PIPE)
        retcode = p.wait()
        # fix known r2r base-pair issue if PDF not created
        if retcode != 0:
            sto = 0
            with open(outfolder + "/" + group + ".sto", 'U') as fin:
                sto = fin.readlines()
            sto[-2] = "# =GF R2R SetDrawingParam autoBreakPairs true\n"
            sto[-1] = "//\n"
            with open("%s/%s.sto" % outinfo, 'w') as fout:
                fout.write(''.join(sto))
            p = Popen(["r2r", "%s/%s.sto" % outinfo, "%s/%s.pdf" % outinfo],
                      stdout=PIPE)
            p.wait()
    except Exception as e:
        print("r2r: ", format_exc(e))


def run_infernal(cmfile, rnd, seqs, outfolder, cpus=1, score=0.0,
                 calibrate=False):
    if exists("%s/R%ihits.fna" % (outfolder, rnd)):
        return
    if not exists(cmfile):
        raise IOError("cmfile path provided does not exist: %s" % cmfile)
    params = {'--mid': True, '--Fmid': 0.0002, '--notrunc': True,
              '--toponly': True, '--cpu': cpus}  # '-g': True,
    if calibrate:
        calibrate_file(cmfile, cpus=cpus)
    result = cmsearch_from_file(cmfile, seqs, RNA, cutoff=score,
                                params=params)
    with open("%s/R%ihits.fna" % (outfolder, rnd), 'w') as fout:
        for hit in result:
            fout.write(">%s score:%0.1f e-val:%f\n%s\n" %
                       (hit[0], hit[14], hit[15], seqs.getSeq(hit[0])))
    if exists("%s/log.txt" % outfolder):
        with open("%s/log.txt" % outfolder, 'a') as fout:
            fout.write("Round %i: %i hits\n" % (rnd, len(result)))


def calculate_overlap(basefolder, rnd):
    basefolder = basefolder.strip()
    if basefolder[:-1] != "/":
        basefolder += "/"

    # get list of groups
    groups = walk(basefolder).next()[1]
    if "fasta_groups" in groups:
        groups.remove("fasta_groups")
    sizegroups = len(groups)

    # compare all group sequences to other sequences
    overlapmatrix = empty(shape=(sizegroups, sizegroups), dtype=float)
    overlapmatrix.fill(-1.0)
    for group, gfolder in enumerate(groups):

        # load in headers of current group
        headers = set([])
        firstgroup = open("%s%s/R%ihits.fna" % (basefolder, gfolder, rnd))
        for header, seq in MinimalFastaParser(firstgroup):
            headers.add(header.split(" ")[0])
        firstgroup.close()
        groupsize = float(len(headers))
        overlapmatrix[group][group] = 1.0
        # go through all other groups and compare sequence headers
        for secgroup in range(group+1, sizegroups):
            secgroupfile = "%s%s/R%ihits.fna" % (basefolder, groups[secgroup],
                                                 rnd)
            if not exists(secgroupfile):
                raise IOError("File not found: %s" % secgroupfile)
            count = 0
            compfile = open(secgroupfile)
            for header, seq in MinimalFastaParser(compfile):
                if header.split(" ")[0] in headers:
                    count += 1
            compfile.close()
            overlapmatrix[group][secgroup] = count / groupsize

    return groups, overlapmatrix


def create_families(basefolder, rnd, outfile=None, fam_cutoff=0.9):
    basefolder = basefolder.strip()
    if basefolder[-1] != "/":
        basefolder += "/"
    groups, overlap = calculate_overlap(basefolder, rnd)
    families = []
    for rowpos, row in enumerate(overlap):
        currgroup = groups[rowpos]
        currpos = -1
        #  find position for family with current group, or new family if needed
        for pos, fam in enumerate(families):
            if currgroup in fam:
                currpos = pos
                break
        if currpos == -1:
            #  not found so new family
            currpos = len(families)
            families.append(set([groups[rowpos]]))
        for colpos, percent in enumerate(row):
            if percent > fam_cutoff:
                families[currpos].add(groups[colpos])
    #  turn sets to lists for returning
    for pos in range(0, len(families)):
        families[pos] = list(families[pos])
    if outfile is not None:
        savetxt(basefolder+outfile, overlap, delimiter="\t",
                header="\t".join(groups))
    return families


def parse_fams_r2r(fam_groups, fam, basefolder, cpus=1):
    """Aligns each family with cm of each group and creates r2r for alignment

    Parameters
    ----------
    fam_groups : list of str
        groups in the family
    fam : str
        name of the family
    basefolder : str
        path to the base output directory
    """
    # Load family stockholm file
    sto = next(StockholmAlignment.from_file(
               join(basefolder, fam, "bayesfold-aln.sto"), RNASequence))
    # grab r2r weights information
    r2r_counts = sto.gf["USE_THIS_WEIGHT_MAP"]
    # write out degapped sequences
    degapped = join(basefolder, fam, "degapped.fna")
    with open(degapped, 'w') as fout:
        fout.write(sto.degapped().to_fasta())
    # apply r2r weights to each family and re-create r2r drawing
    for group in fam_groups:
        # align family sequences to cm for group
        group_sto = join(basefolder, fam, "%s.sto" % group)
        cmalign(degapped, join(basefolder, group, "cmfile.cm"),
                group_sto, cpus)
        # add weight to stockholm file and write back out
        sto = StockholmAlignment.from_file(group_sto, RNASequence)
        sto.gf["USE_THIS_WEIGHT_MAP"] = r2r_counts
        with open(group_sto, 'w') as fout:
            fout.write(str(sto))
        # rebuild r2r with new alignment
        make_r2r(group_sto, join(basefolder, fam), "%s_%s" % (fam, group))
