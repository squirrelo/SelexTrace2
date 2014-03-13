from sys import stdout
from os.path import exists, dirname, abspath
from os import mkdir
from subprocess import Popen, PIPE
from math import ceil
from random import shuffle
from multiprocessing import Pool, Manager
from traceback import format_exc

from cogent import LoadSeqs, RNA
from cogent.core.sequence import RnaSequence
from cogent.app.infernal_v11 import (cmsearch_from_file, cmbuild_from_file,
                                     calibrate_file)
from cogent.parse.fasta import MinimalFastaParser
from cogent.format.stockholm import stockholm_from_alignment
from nwalign import global_align, score_alignment

from bayeswrapper import bayesfold
from selextrace.stutils import count_seqs


def fold_clusters(lock, cluster, seqs, otufolder):
    '''Function for multithreading.
    Computes structure for a cluster and writes it to file'''
    aln, struct = bayesfold(seqs, params={"-diags": True})
    #write structure out to file
    try:
        lock.acquire()
        cfo = open(otufolder + "cluster_structs.fasta", 'a')
        cfo.write(">%s\n%s\n" % (cluster, struct))
        cfo.write(aln.toFasta() + "\n")
        cfo.close()
        lock.release()
    except Exception:
        lock.release()


class SeqStructure(object):
    def __init__(self, structure, seq=None):
        self.struct = structure
        self.seq = None
        self._structmap = self._create_structmap(structure)
        self._seqmap = None
        if seq is not None:
            self.add_seqmap(seq)
            self.seq = seq

        self.pairs = ('AU', 'UA', 'GC', 'CG', 'GU', 'UG')

    def __len__(self):
        return len(self.struct)

    def __str__(self):
        if self._seqmap is not None:
            strmap = "seqmap:" + ' '.join([str(pos) for pos in self._seqmap])
        else:
            strmap = "structmap: " + ' '.join([str(pos)
                                               for pos in self._structmap])
        return ''.join([self.struct, "\n", strmap])

    def _create_structmap(self, structure):
        stack = []
        structmap = []
        for pos, symbol in enumerate(structure):
            if symbol == "(":
                stack.append(pos)
            elif symbol == ")":
                p1 = stack.pop()
                structmap.append((p1, pos))
        return structmap

    def add_seqmap(self, seq):
        if len(seq) != len(self.struct):
            raise ValueError("sequence must be same length as structure!")
        seqmap = []
        for p1, p2 in self._structmap:
            pair = seq[p1] + seq[p2]
            seqmap.append((p1, p2, pair))
        self._seqmap = seqmap
        self.seq = seq

    def score_seq(self, seq):
        """Scores sequence based on how well it fits into given structure.
            INPUT:
            seq: RNA sequence to thread into structure
            OUTPUT
            Score, where pair is given 2 points if matches pair in structure
            and one point if it can base pair at all. This divided by seq
            length to get final score between 0 and 1.
        """
        maxscore = 0
        score = 0
        start = 0
        end = len(seq)
        if self._seqmap is None:
            raise RuntimeError("No seqmap exists!")
        if len(seq) == len(self.struct):
        #same length so trivial
            maxscore = self._eval_struct_seq(seq)
        elif len(seq) > len(self.struct):
            #sequence longer than struct so keep slicing sequence
            #evaluate all possible sequence fittings into the structure
            #return highest scoring one
            while end != len(seq):
                score = self._eval_struct_seq(seq[start:end])
                if maxscore < score:
                    maxscore = score
                start += 1
                end += 1
        else:
            #struct longer than sequence so walk seq down struct and compare
            #evaluate all possible sequence fittings into the structure
            #return highest scoring one
            while end != len(self.struct):
                score = self._eval_struct_seq(seq, seqoffset=start)
                if maxscore < score:
                    maxscore = score
                start += 1
                end += 1
        #return normalized score
        return float(maxscore) / len(self._seqmap)

    def _eval_struct_seq(self, seq, seqoffset=0):
        if self._seqmap is None:
            raise RuntimeError("There is no seqmap to evaluate!")
        score = 0
        for p1, p2, pair in self._seqmap:
            p2 = p2 - seqoffset
            #make sure we are in range of sequence given including offset
            if p1 < seqoffset or p2 >= len(seq):
                continue
            p1 = p1 - seqoffset
            #p1 & p2: positions of bases, pair: basepair at that position
            currpair = seq[p1] + seq[p2]
            currpair = currpair.upper()
            pair = pair.upper()
            if currpair in self.pairs:
                if currpair == pair:
                    #perfect pair match with bases scores 2 points
                    score += 2
                else:
                    #ability to basepair at all scores 1 point
                    score += 1
        return score


def build_reference(keys, refsize):
    '''Creates a random list of references refsize long, returning rest as
       nonref
    '''
    shuffle(keys)
    #return reference, nonreference by slicing list
    return keys[:refsize], keys[refsize:]


def group_to_reference(reference, nonref, minscore, cpus=1):
    pool = Pool(processes=cpus)
    manager = Manager()
    groupstruct = manager.dict()
    nogroup = manager.list()

    chunksize = len(nonref)/cpus
    if chunksize == 0:
        chunksize = 1
    for startpos in range(0, len(nonref), chunksize):
        #divide up nonref into chunks and align each chunk to reference seqs
        #final # chunks == number of cpus available
        endpos = startpos+chunksize
        if endpos >= len(nonref):
            pool.apply_async(func=group, args=(nonref[startpos:], minscore,
                                               reference, groupstruct,
                                               nogroup))
            break
        pool.apply_async(func=group, args=(nonref[startpos:endpos],
                                           minscore, reference, groupstruct,
                                           nogroup))
    pool.close()
    pool.join()
    #need to explicitly change to list and dict because manager versions are
    #missing certain key functionality expected of lists and dicts
    return dict(groupstruct), list(nogroup)


def group(nonref, minscore, ref=None, groupstruct=None, nogroup=None):
    #takes in list of seqstructure objects for nonref and ref
    try:
        basefolder = dirname(abspath(__file__))
        matrix = basefolder + "/NucMatrix"
        denovo = True
        #if ref list is pased, know we are refrence grouping
        if ref is not None:
            denovo = False
        if groupstruct is None:
            groupstruct = {}
        if nogroup is None:
            nogroup = []
        #loop through all nonreference items
        for pos, currnonref in enumerate(nonref):
            seq1 = currnonref.seq.replace("-", "")
            bestref = ""
            bestscore = minscore
            if denovo:
                ref = nonref[pos+1:]
            #compare to each reference item
            for teststruct in ref:
                seq2 = teststruct.seq.replace("-", "")
                #get alignment score and add to seq/struct score
                aln = global_align(seq1, seq2, gap_open=-1, gap_extend=-1,
                                   matrix="selextrace/NucMatrix")
                alnsc = score_alignment(aln[0], aln[1], gap_open=-1, 
                                        gap_extend=-1, matrix=matrix)
                #score is normalized by dividing each score by sequence length
                #then adding. This should keep scores between zero and one
                score = (alnsc/len(aln[0]) + currnonref.score_seq(seq2))/3
                if score >= bestscore:
                    bestscore = score
                    bestref = teststruct.struct
            if bestref != "":
                if bestref not in groupstruct:
                    groupstruct[bestref] = [currnonref.struct]
                else:
                    groupstruct[bestref].append(currnonref.struct)
            else:
                nogroup.append(currnonref)
    except Exception, e:
        print "GROUP: ", format_exc(e)
    return groupstruct, nogroup


def group_by_seqstruct(structgroups, structscore, specstructs=None,
                       setpercent=0.01, cpus=1):
        '''Does grouping by way of de-novo reference creation and clustering
            structgroups - dictionary with ALL structures and the Alignment
                           object keyed to them
            structscore - maximum score to consider grouping structures
            specstructs - a list of a subset of structures in structgroups
                          to cluster (optional)
            setpercent - Allows manual percentage setting for ref structures
                           (default 1%  of dict)
        '''
        #fail if nothing to compare
        if len(structgroups) < 1:
            raise ValueError("Must have at least one structure to group!")
        #return the list directly if only one item (useful for breakout work)
        if len(structgroups) == 1:
            return {structgroups.keys()[0]: []}
        grouping = []
        if specstructs is None:
            for currstruct in structgroups:
                seq = ''.join(structgroups[currstruct].majorityConsensus())
                grouping.append(SeqStructure(currstruct, seq))
        else:
            for currstruct in specstructs:
                seq = ''.join(structgroups[currstruct].majorityConsensus())
                grouping.append(SeqStructure(currstruct, seq))
        #create SeqStructures objects for items we are clustering
        #just de-novo group if 20 or less to save time and effort
        if len(grouping) <= 20:
            grouped, ungrouped = group(grouping, structscore)
            for ug in ungrouped:
                grouped[ug.struct] = []
            return grouped
        #for speed, get 1% as initial clustering or user defined.
        #Need at least 10 structs though.
        finishlen = int(ceil(len(grouping) * setpercent))
        if finishlen < 10:
            finishlen = 10
        startgrouped = 0
        endgrouped = 1
        grouped = {}
        ungrouped = grouping
        # keep refining while not at limit and are still grouping structs
        while (len(grouping) - len(grouped) > finishlen and
               startgrouped != endgrouped):
            startgrouped = len(grouped)
            ref, nonref = build_reference(ungrouped, finishlen)
            g, ungrouped = group_to_reference(ref, nonref, structscore, cpus)
            grouped.update(g)
            endgrouped = len(grouped)
        #do the last grouping
        g, ungrouped = group(ungrouped, structscore)
        grouped.update(g)
        #add ungroupable bit to end
        for ug in ungrouped:
            grouped[ug.struct] = []
        return grouped


def final_fold(seqs, params, outfolder, group, fold=True):
            if exists("%sgroup_%i.fasta" % (outfolder, group)):
                return
            if fold:
                aln, struct = bayesfold(seqs, params=params)
            else:
                aln = LoadSeqs(data=seqs, moltype=RNA)
            aln.Names.sort(reverse=True, key=lambda c: count_seqs(c))
            fout = open("%sgroup_%i.fasta" % (outfolder, group), 'w')
            fout.write(aln.toFasta() + "\n")
            if fold:
                fout.write(">SS_Struct\n%s\n" % struct)
            fout.close()


def create_group_output(groupfasta, basefolder, minseqs=1, cpus=1):
    '''Function for multithreading. Creates the final BayesFold alignment and
    writes to files, then r2r struct and infernal CM file'''
    try:
        #skip if already run and program just crashed or whatever
        currgroup = groupfasta.split("/")[-1].split(".")[0]
        currotufolder = basefolder + currgroup
        if exists(currotufolder):
            return
        seqs = []
        weights = []
        maxweight = 0
        count = 0
        fin = open(groupfasta, 'rU')
        for header, seq in MinimalFastaParser(fin):
            if header == "SS_Struct":
                struct = seq
                continue
            seqs.append((header.split()[0], seq.strip()))
            weight = count_seqs(header)
            count += weight
            if weight > maxweight:
                maxweight = weight
            weights.append(header.split()[0])
            weights.append(str(weight))
        fin.close()
        aln = LoadSeqs(data=seqs, moltype=RNA)
        del seqs
        if count < minseqs:
            return
        mkdir(currotufolder)
        out = ' '.join([currgroup, ":\n", str(count),
                        "sequences\n", str(aln.getNumSeqs()),
                        "unique sequences\nStructure: ", struct, "\n"])
        #write out alignment and structure in fasta format
        logout = open(currotufolder + "/log.txt", 'w')
        logout.write(out)
        logout.close()
        with open(currotufolder + "/bayesfold-aln.fasta", 'w') as alnout:
            alnout.write(">SS_struct\n%s\n%s" % (struct, aln.toFasta()))

        #create standard weights for infernal
        infweights = ""
        for pos in range(0, len(weights), 2):
            infweights = ''.join([infweights, '#GS\t%s\tWT\t%s\n' %
                                 (weights[pos],
                                  str(float(weights[pos+1]) / maxweight))])
        #create weights in for r2r
        r2r_weights = "#=GF USE THIS WEIGHT MAP " + ' '.join(weights)
        #create sto file with r2r and std weights
        sto = stockholm_from_alignment(aln, GC_annotation={'SS_cons': struct})
        sto = sto.split("\n")
        sto[-1] = infweights
        sto.append(r2r_weights + "\n")
        sto.append("//\n")
        stofile = currotufolder + "/bayesfold-aln.sto"
        with open(stofile, 'w') as alnout:
            alnout.write('\n'.join(sto))

        #make R2R secondary structure for alignment
        make_r2r(currotufolder+"/bayesfold-aln.sto", currotufolder, currgroup)
        #create CM file for infernal from group
        with open(currotufolder + "/cmfile.cm", 'w') as fout:
            fout.write(cmbuild_from_file(stofile, params={'--wgiven': True}))
        calibrate_file(currotufolder + "cmfile.cm", cpus=cpus)
    except Exception, e:
        print "create_group_output:\n", format_exc(e)
        stdout.flush()


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
        #fix known r2r base-pair issue if PDF not created
        if retcode != 0:
            sto = 0
            with open(outfolder + "/" + group + ".sto", 'U') as fin:
                sto = fin.readlines()
            sto[-2] = "#=GF R2R SetDrawingParam autoBreakPairs true\n"
            sto[-1] = "//\n"
            with open("%s/%s.sto" % outinfo, 'w') as fout:
                fout.write(''.join(sto))
            p = Popen(["r2r", "%s/%s.sto" % outinfo, "%s/%s.pdf" % outinfo],
                      stdout=PIPE)
            p.wait()
    except Exception, e:
        print "r2r: ", format_exc(e)


def score_local_rnaforester(struct1, struct2):
    '''returns local aignment score of two structures'''
    #return gigantically negative number if no structure for one struct
    if "(" not in struct1 or "(" not in struct2:
        raise ValueError("%s\n%s\nNo pairing in structures!" % (struct1,
                         struct2))
    p = Popen(["RNAforester", "--score", "-l"], stdin=PIPE, stdout=PIPE)
    p.stdin.write(''.join([struct1, "\n", struct2, "\n&"]))
    return int(p.communicate()[0].split("\n")[-2])


def score_multi_forester(basestruct, checkstruct, foresterscore):
    if score_local_rnaforester(basestruct, checkstruct) > foresterscore:
        return checkstruct
    else:
        return ""


def group_by_forester(fulldict, foresterscore, cpus=1):
    structs = fulldict.keys()
    for pos, currstruct in enumerate(structs):  # for each structure
        #skip if already grouped
        if currstruct not in fulldict:
            continue
        scores = set([])
        pool = Pool(processes=cpus)
        #compare everything as fast as possible using multiprocessing
        #comparisons end up as set of structs above threshold plus ""
        for teststruct in structs[pos+1:]:
            pool.apply_async(func=score_multi_forester,
                             args=(currstruct, teststruct, foresterscore),
                             callback=scores.add)
        pool.close()
        pool.join()
        #remove empty and add remaining structs to currgroup
        if "" in scores:
            scores.discard("")
        for struct in scores:
            if struct in fulldict:
                fulldict[currstruct].extend(fulldict[struct])
                fulldict.pop(struct)
    return fulldict


def run_infernal(cmfile, rnd, basefolder, outfolder, cpus=1, score=0.0):
    seqs = 0
    #Only search unique sequences to save time
    #check if previous run has removed some sequences, load correct file
    if not exists(cmfile):
        raise IOError("cmfile path provided does not exist: %s" % cmfile)
    uniques_remain_file = ''.join([basefolder, "R", str(rnd), "/R", str(rnd),
                                   "-Unique-Remaining.fasta"])
    uniques_file = "%sR%i/R%i-Unique.fasta" % (basefolder, rnd, rnd)
    if exists(uniques_remain_file):
        seqs = LoadSeqs(uniques_remain_file, moltype=RNA, aligned=False)
    elif exists(uniques_file):
        seqs = LoadSeqs(uniques_file, moltype=RNA, aligned=False)
    else:
        raise IOError("Round's fasta file does not exist!")
    params = {'--mid': True, '--Fmid': 0.0002, '--notrunc': True,
              '--toponly': True, '--cpu': cpus}  # '-g': True,
    result = cmsearch_from_file(cmfile, seqs, RNA, cutoff=score,
                                params=params)
    with open("%s/R%ihits.fna" % (outfolder, rnd), 'w') as fout:
        for hit in result:
            fout.write(">%s score:%0.1f e-val:%f\n%s\n" % (hit[0], hit[14],
                                                           hit[15],
                                                           seqs.getSeq(hit[0])))
    if exists("%s/log.txt" % outfolder):
        with open("%s/log.txt" % outfolder, 'a') as fout:
            fout.write("Round %i: %i hits\n" % (rnd, len(result)))
