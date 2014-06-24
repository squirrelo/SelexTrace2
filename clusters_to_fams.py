# !/usr/bin/env python

__author__ = "Joshua Shorenstein"
__copyright__ = "Copyright 2014, SelexTrace project"
__credits__ = ["Joshua Shorenstein"]
__license__ = "BSD"
__version__ = "0.0.1-dev"
__maintainer__ = "Joshua Shorenstein"
__email__ = "joshua.shorenstein@colorado.edu"
__status__ = "Development"

from os.path import exists, join
from os import mkdir, walk
import argparse
from time import time
from datetime import datetime
from multiprocessing import Pool, Manager

from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs, RNA

from selextrace.stutils import cluster_seqs, count_seqs
from selextrace.ctilib import (fold_clusters, create_final_output,
                               group_by_seqstruct, align_order_seqs,
                               run_infernal, create_families,
                               write_clusters, read_clusters,
                               create_seqstructs, parse_fams_r2r)

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
    parser.add_argument('--minseqs', type=int, default=-1,
        help=("Min number of seqs for group to be significant "
              "(Default 0.1% total)"))
    parser.add_argument('--csc', type=float, default=0.8,
        help="Score cutoff for clustering (Default 0.8)")
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

    #calculate minseqs if necessary
    if args.minseqs == -1:
        with open(args.i) as fin:
            args.minseqs = int(count_seqs([h for h, s in
                                           MinimalFastaParser(fin)]) * 0.001)

    if basefolder[:-1] != "/":
        basefolder += "/"
    if outfolder[:-1] != "/":
        outfolder += "/"
    if not exists(outfolder):
        mkdir(outfolder)

    date = str(datetime.now())
    print "Program started ", date

    # print out run info to a file
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
    secs = time()
    clustfile = outfolder + "clusters.txt"
    structfile = outfolder + "cluster_structs.fasta"

    if exists(structfile):
        # dont need to do anything since already folded by next step
        print "Sequences previously clustered"
        with open(clustfile) as fin:
            numclusts = int(fin.readline().strip())
    elif exists(clustfile):
        # already clustered but not folded, so read in clusters
        with open(clustfile) as fin:
            clusters, numclusts = read_clusters(fin)

        print "Sequences previously clustered, %i clusters" % numclusts
    else:
        print "Running uclust over sequences"
        # cluster the initial sequences by sequence simmilarity
        clusters = cluster_seqs(args.i, args.sim, folderout=args.o,
                                gapopen='1.0', gapext='1.0')
        with open(clustfile, 'w') as fout:
            write_clusters(clusters, fout)
        numclusts = len(clusters)
        clusters.clear()
        del clusters
        print "Runtime: %0.2f min" % ((time() - secs)/60)

    if not exists(structfile):
        # create file to write to if not already there
        with open(structfile, 'w'):
            pass
        print "Running BayesFold over %i clusters" % numclusts
        secs = time()
        # make a pool of workers, one for each cpu available
        manager = Manager()
        pool = Pool(processes=args.c)
        lock = manager.Lock()
        with open(clustfile) as fin:
            # run the pool over all clusters to get file of structures
            # throw out first line since it has number in it
            fin.readline()
            currclust = fin.readline().strip(">").strip()
            fin.readline()
            cluster = []
            for header, seq in MinimalFastaParser(fin):
                if header.startswith("cluster_"):
                    pool.apply_async(func=fold_clusters, args=(lock, currclust,
                                     cluster, structfile))
                    currclust = header
                    cluster = []
                else:
                    cluster.append((header, seq))
            pool.apply_async(func=fold_clusters, args=(lock, currclust,
                             cluster, structfile))
        pool.close()
        pool.join()
        cluster = None
        del cluster
    else:
        print "Clusters previously folded"

    # read in all structures now that they are folded and aligned
    with open(structfile) as fin:
        seqstructs = create_seqstructs(fin, numclusts)

    print "Runtime: %0.2f min" % ((time() - secs)/60)
    print "==Grouping clusters by sequence & secondary structure=="
    if not exists(outfolder + "fasta_groups/"):
        secs = time()
        print "start: %i initial groups" % len(seqstructs)
        # initial clustering by structures generated in first folding
        # run the pool over all shape groups to get final grouped structgroups
        grouped = group_by_seqstruct(seqstructs, clustscore, cpus=args.c,
                                     setpercent=0.01)
        del seqstructs

        print "%i end groups (%0.2f hrs)" % (len(grouped), (time()-secs)/3600)
        print "Align end groups"

        secs = time()
        # collect all structs together and post-process them.
        # Does final fold and file printout
        with open(clustfile) as fin:
            clusters, numclusts = read_clusters(fin)
        mkdir(outfolder + "fasta_groups")
        params = {"-diags": True, "-maxiters": 5}
        pool = Pool(processes=args.c)
        # need to use fasta string because pycogent sequence collections
        # HATE multithreading so can't use them
        for num, ref in enumerate(grouped):
            seqs = clusters[ref]
            fastafolder = outfolder+"/fasta_groups/"
            for cluster in grouped[ref]:
                seqs.extend(clusters[cluster])
            pool.apply_async(func=align_order_seqs, args=(seqs, params,
                                                          fastafolder, num))
        grouped.clear()
        del grouped
        clusters.clear()
        del clusters
        pool.close()
        pool.join()
        print "Align complete (%0.2f min)" % ((time()-secs)/60)
    else:
        print "Previously grouped"

    print "==Creating final groups=="
    secs = time()
    # smaller pool for memeory savings, 1 task per child to try and gc each run
    # NEED TO FIX BAYESFOLD ARRAY2D BEING A MEMORY HOG
    procs = args.c / 4
    if procs < 1:
        procs = 1
    cpus = args.c / procs
    if cpus < 1:
        cpus = 1
    pool = Pool(processes=procs, maxtasksperchild=1)
    # run the pool over all groups to get final structures
    for group in walk(outfolder + "fasta_groups").next()[2]:
        pool.apply_async(func=create_final_output,
                         args=(outfolder+"fasta_groups/"+group, outfolder,
                               args.minseqs, cpus))
    pool.close()
    pool.join()

    # get sequence counts for each group
    grouporder = []
    groups = {}
    count = 0
    groupslist = [g for g in walk(outfolder).next()[1] if "group_" in g]
    for group in groupslist:
        if group == "fasta_groups":
            continue
        count += 1
        # read in group sequence counts
        log = open(outfolder + group + "/log.txt")
        loginfo = log.readlines()
        log.close()
        grouporder.append((group, int(loginfo[1].split()[0]),
                          int(loginfo[2].split()[0])))
        # read in group structure and build dict for families creation
        with open(outfolder + group + "/bayesfold-aln.fasta") as structin:
            struct = structin.readlines()[-2].strip()
        groups[struct] = [group]

    # write out file of sequence counts
    grouporder.sort(reverse=True, key=lambda x: x[1])
    groupsizefile = open(outfolder + "/group_sizes.txt", 'w')
    groupsizefile.write("Group\tTotal Seqs\tUnique Seqs\n")
    for info in grouporder:
        groupsizefile.write("%s\t%s\t%s\n" % info)
    groupsizefile.close()

    print count, "final groups"
    print "Runtime: %0.2f hrs" % ((time() - secs) / 3600)

    print "==Creating families from groups=="
    secs = time()
    groups = walk(outfolder).next()[1]
    uniques_remain_file = "%sR%i/R%i-Unique-Remaining.fasta" % (basefolder,
                                                                args.r, args.r)
    uniques_file = "%sR%i/R%i-Unique.fasta" % (basefolder, args.r, args.r)
    if exists(uniques_remain_file):
        seqs = LoadSeqs(uniques_remain_file, moltype=RNA, aligned=False)
    elif exists(uniques_file):
        seqs = LoadSeqs(uniques_file, moltype=RNA, aligned=False)
    else:
        raise IOError("Round's fasta file does not exist!")
    for group in groups:
        # run infernal over current round for all groups
        if group == "fasta_groups":
            continue
        run_infernal("%s%s/cmfile.cm" % (outfolder, group), args.r, seqs,
                     outfolder + group + "/", cpus=args.c, score=args.isc)
    del seqs

    # get families and write out families list
    families = create_families(outfolder, args.r, outfile="family_matrix.txt")
    print len(families), "families"
    with open(outfolder + "families.txt", 'w') as fout:
        for famnum, family in enumerate(families):
            fout.write("fam_%i\n" % famnum)
            fout.write('\n'.join(family) + '\n')

    fambase = outfolder + "fasta_families/"
    mkdir(fambase)
    for famnum, family in enumerate(families):
        # combine all seqs for family. Abusing alignment obj for removing gaps
        # need to use fasta because pycogent addSeqs does not cooperate if not
        # alignment and all seqs same length
        famseqs = LoadSeqs("%sfasta_groups/%s.fna" % (outfolder, family[0]),
                           moltype=RNA).degap().toFasta() + "\n"
        for group in family[1:]:
            groupseqsfile = "%sfasta_groups/%s.fna" % (outfolder, group)
            if not exists(groupseqsfile):
                raise IOError("FASTA NOT FOUND: %s" % groupseqsfile)
            famseqs += LoadSeqs(groupseqsfile, moltype=RNA).degap().toFasta()\
                + "\n"
        # fold family sequences
        params = {"-diags": True, "-maxiters": 5}
        align_order_seqs(famseqs, params, fambase, famnum, prefix="fam_")
        create_final_output("%sfam_%i.fna" % (fambase, famnum), outfolder,
                            args.minseqs, args.c)

        # use r2r output and map to each group in family
        with open(join(outfolder, "families.txt")) as fin:
            fam_groups = []
            currfam = fin.readline().strip()
            for line in fin:
                line = line.strip()
                if "fam_" in line:
                    parse_fams_r2r(fam_groups, currfam, cpus=args.c)
                    currfam = line
                    fam_groups = []
                else:
                    fam_groups.append(line)

    print "Runtime: %0.2f hrs" % ((time() - secs) / 3600)

    endtime = (time() - starttime)/3600
    print "Program ended", datetime.now(), " Runtime:", endtime, "hrs"
