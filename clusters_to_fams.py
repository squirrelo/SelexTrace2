#!/usr/bin/env python

__author__ = "Joshua Shorenstein"
__copyright__ = "Copyright 2014, SelexTrace project"
__credits__ = ["Joshua Shorenstein"]
__license__ = "BSD"
__version__ = "0.0.1-dev"
__maintainer__ = "Joshua Shorenstein"
__email__ = "joshua.shorenstein@colorado.edu"
__status__ = "Development"

from os.path import exists
from os import mkdir, walk
import argparse
from time import time
from datetime import datetime
from multiprocessing import Pool, Manager

from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs, RNA

from selextrace.stutils import cluster_seqs, count_seqs
from selextrace.ctilib import (fold_clusters, create_group_output,
                               group_by_seqstruct, group_by_forester,
                               align_order_seqs, run_infernal)

if __name__ == "__main__":
    starttime = time()
    parser = argparse.ArgumentParser(description="Runs sequence grouping \
        and infernal over all rounds of a SELEX selection")
    parser.add_argument('-i', required=True,
        help="FASTA file of unique sequences sorted by abundance.")
    parser.add_argument('-o', required=True,
        help="Base folder to output all data")
    parser.add_argument('-r', required=True, type=int,
        help="Round this input file comes from")
    parser.add_argument('-f', required=True,
        help="Base folder containing SELEX round fasta files")
    parser.add_argument('--sim', type=float, default=0.99,
        help="Simmilarity for uclust. (Default 0.99)")
    parser.add_argument('--minseqs', type=int, default=100,
        help="Min number of seqs for group to be significant (Default 100)")
    parser.add_argument('--csc', type=float, default=0.75,
        help="Score cutoff for clustering (Default 0.75)")
    parser.add_argument('--isc', type=int, default=80,
        help="Score cutoff for Infernal. (Default 80)")
    parser.add_argument('--fsc', type=int, default=150,
        help="Score cutoff for RNAforester. (Default 150)")
    parser.add_argument('-c', type=int, default=1,
        help="Number of CPUs to use (Default 1)")

    args = parser.parse_args()
    if args.c < 1:
        raise ValueError("ERROR: CPU count must be at least 1!")
    if args.fsc < 0:
        raise ValueError("ERROR: RNAforester score cutoff must be >0!")
    if args.sim <= 0.0 or args.sim > 1.0:
        raise ValueError("ERROR: clustering simmilarity must be > 0 and <= 1!")
    clustscore = args.csc
    outfolder = args.o.strip()
    basefolder = args.f.strip()

    if not exists(basefolder):
        raise IOError("Basefolder does not exist!")

    if basefolder[:-1] != "/":
        basefolder += "/"
    if outfolder[:-1] != "/":
        outfolder += "/"
    if not exists(outfolder):
        mkdir(outfolder)

    date = str(datetime.now())
    print "Program started ", date

    #print out run info to a file
    infofile = open(outfolder + "runparams.txt", 'w')
    infofile.write(''.join(["Program started ", date, "\n",
                            "FASTA file:\t", args.i, "\n",
                            "Output folder:\t", args.o, "\n"
                            "Round:\t", str(args.r), "\n"
                            "Uclust simmilarity:\t", str(args.sim), "\n",
                            "Min seqs for group:\t", str(args.minseqs), "\n",
                            "Clustering score cutoff:\t", str(args.csc), "\n",
                            "Infernal score cutoff:\t", str(args.isc), "\n",
                            "RNAforester score cutoff:\t", str(args.fsc), "\n",
                            "CPUs:\t", str(args.c), "\n"]))
    infofile.close()
    print "==Clustering sequences by primary sequence=="
    clusters = {}
    secs = time()
    if exists(outfolder + "cluster_structs.fasta"):
        #dont need to do anything since already folded by next step
        print "Sequences previously clustered"
        cin = open(outfolder + "clusters.txt")
        numclusts = int(cin.readline().strip())
        cin.close()
    elif exists(outfolder + "clusters.txt"):
        #already clustered but not folded, so read in clusters
        clustersin = open(outfolder + "clusters.txt")
        numclusts = int(clustersin.readline().strip())
        currclust = ""
        for header, seq in MinimalFastaParser(clustersin):
            if "cluster_" in header:
                currclust = header
                clusters[currclust] = []
            else:
                clusters[currclust].append((header, seq))
        clustersin.close()
        print "Sequences previously clustered,", numclusts, "clusters"
    else:
        print "Running uclust over sequences"
        #cluster the initial sequences by sequence simmilarity
        clusters = cluster_seqs(args.i, args.sim, folderout=args.o,
                                gapopen='1.0', gapext='1.0')

        #print that shit to file
        hold = clusters.keys()
        hold.sort()
        cout = open(outfolder + "clusters.txt", 'w')
        cout.write(str(len(clusters)) + "\n")
        for cluster in hold:
            cout.write(">%s\n%s\n" % (cluster, cluster))
            for header, seq in clusters[cluster]:
                cout.write(">%s\n%s\n" % (header, seq))
        cout.close()
        numclusts = len(clusters)
        print str(len(clusters)) + " clusters"
        print "Runtime: %0.2f min" % ((time() - secs)/60)

    if not exists(outfolder + "cluster_structs.fasta"):
        #create file to write to if not already there
        with open(outfolder + "cluster_structs.fasta", 'w'):
            pass
        print "Running BayesFold over " + str(len(clusters)) + " clusters"
        secs = time()
        #make a pool of workers, one for each cpu available
        #200 tasks per child to force garbage collection, since OSX sucks at it
        manager = Manager()
        pool = Pool(processes=args.c, maxtasksperchild=200)
        lock = manager.Lock()
        #run the pool over all clusters to get file of structures
        for cluster in clusters:
            pool.apply_async(func=fold_clusters, args=(lock, cluster,
                             clusters[cluster], outfolder))
        pool.close()
        pool.join()
    else:
        print "Clusters previously folded"

     #read in all structures now that they are folded and aligned
    structgroups = {}
    count = 1
    cfo = open(outfolder + "cluster_structs.fasta", 'rU')
    #throw out first line (dont need cluster number), then read in currstruct
    cfo.readline()
    currstruct = cfo.readline().strip()
    currseqs = []
    for header, seq in MinimalFastaParser(cfo):
        if "cluster_" in header:
            count += 1
            if currstruct not in structgroups:
            #turn list of tuples into alignment object
                structgroups[currstruct] = LoadSeqs(data=currseqs, moltype=RNA)
            else:
                aln = LoadSeqs(data=currseqs, moltype=RNA)
                structgroups[currstruct] = structgroups[currstruct].addSeqs(aln)
            #move on to next structgroup
            currstruct = seq
            currseqs = []
        else:
            currseqs.append((header, seq))
    structgroups[currstruct] = LoadSeqs(data=currseqs, moltype=RNA)
    cfo.close()
    del currseqs
    for struct in structgroups:
        if len(struct) != len(structgroups[struct]):
            print struct

    if count != numclusts:
        raise AssertionError(" ".join([str(count), "structures,",
                                       str(len(clusters)),
                                       "clusters. Not all clusters folded!"]))
    clusters.clear()
    del clusters

    print "Runtime: %0.2f min" % ((time() - secs)/60)
    print "==Grouping clusters by sequence & secondary structure=="
    if not exists(outfolder + "fasta_groups/"):
        secs = time()
        print "start: " + str(len(structgroups)) + " initial groups"
        #initial clustering by structures generated in first folding
        #run the pool over all shape groups to get final grouped structgroups
        hold = group_by_seqstruct(structgroups, clustscore, cpus=args.c,
                                  setpercent=0.01)
        print "%i end groups (%0.2f hrs)" % (len(hold), (time()-secs)/3600)
        print "Align end groups"

        secs = time()
        #collect all structs together and post-process them.
        #Does final fold and file printout
        mkdir(outfolder + "fasta_groups")
        params = {"-diags": True, "-maxiters": 5}
        pool = Pool(processes=args.c)
        #need to use fasta string because pycogent sequence collections
        #HATE multithreading so can't use them
        for num, struct in enumerate(hold):
            count = 0
            seqs = structgroups[struct].degap().toFasta() + "\n"
            count += count_seqs(structgroups[struct].Names)
            fastafolder = outfolder+"/fasta_groups/"
            for substruct in hold[struct]:
                seqs = ''.join([seqs,
                                structgroups[substruct].degap().toFasta(),
                                "\n"])
                count += count_seqs(structgroups[substruct].Names)
            pool.apply_async(func=align_order_seqs, args=(seqs, params,
                                                          fastafolder, num))
            #final_fold(seqs, params, outfolder + "fasta_groups/", num)
        hold.clear()
        del hold
        pool.close()
        pool.join()
        print "Align complete (%0.2f min)" % ((time()-secs)/60)
    else:
        print "Previously grouped"

    print "==Creating final groups=="
    secs = time()
    #smaller pool for memeory savings, 1 task per child to try and gc each run
    #NEED TO FIX BAYESFOLD ARRAY2D BEING A MEMORY HOG
    procs = args.c / 4
    if procs < 1:
        procs = 1
    cpus = args.c / procs
    if cpus < 1:
        cpus = 1
    pool = Pool(processes=procs, maxtasksperchild=1)
    #run the pool over all groups to get final structures
    for group in walk(outfolder + "fasta_groups").next()[2]:
        pool.apply_async(func=create_group_output,
                         args=(outfolder+"fasta_groups/"+group, outfolder,
                               args.minseqs, cpus))
    pool.close()
    pool.join()

    #get sequence counts for each group
    grouporder = []
    groups = {}
    count = 0
    for group in walk(outfolder).next()[1]:
        if group == "fasta_groups":
            continue
        count += 1
        #read in group sequence counts
        log = open(outfolder + group + "/log.txt")
        loginfo = log.readlines()
        log.close()
        grouporder.append((group, int(loginfo[1].split()[0]),
                          int(loginfo[2].split()[0])))
        #read in group structure and build dict for families creation
        with open(outfolder + group + "/bayesfold-aln.fasta") as structin:
            struct = structin.readlines()[-2].strip()
        groups[struct] = [group]

    #write out file of sequence counts
    grouporder.sort(reverse=True, key=lambda x: x[1])
    groupsizefile = open(outfolder + "/group_sizes.txt", 'w')
    groupsizefile.write("Group\tTotal Seqs\tUnique Seqs\n")
    for info in grouporder:
        groupsizefile.write("%s\t%s\t%s\n" % info)
    groupsizefile.close()

    print count, "final groups"
    print "Runtime: %0.2f hrs" % ((time() - secs) / 3600)

    print "===Running Infernal for groups==="
    print "infernal score cutoff: ", args.isc

    for group in walk(outfolder).next()[1]:
        if group == "fasta_groups":
            continue
        for rnd in range(1, args.r+1):
            run_infernal("%s%s/cmfile.cm" % (outfolder, group), rnd,
                         basefolder, outfolder + group+"/", cpus=args.c,
                         score=args.isc)

    print "==Creating families from groups=="
    print "RNAforester score cutoff:", args.fsc
    print len(groups), "initial groups"
    #Now group by forester local aignment to get larger families grouped
    groups = group_by_forester(groups, args.fsc, args.c)
    fout = open(outfolder + "families.txt", 'w')
    for fam in groups:
        fout.write("\t".join(groups[fam]) + "\n")
    fout.close()
    print len(groups), "final families"
    print "Runtime: " + str((time() - secs) / 3600) + " hrs"
    groups.clear()
    del groups

    endtime = (time() - starttime)/3600
    print "Program ended", datetime.now(), " Runtime:", endtime, "hrs"
