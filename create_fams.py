import argparse
from time import time
from os import mkdir, walk
from os.path import exists

from cogent import LoadSeqs, RNA

from selextrace.ctilib import (create_families, run_infernal, align_order_seqs,
                               create_final_output)

if __name__ == "__main__":
    secs = time()
    parser = argparse.ArgumentParser(description="Runs family creation on"
        "previously grouped SelexTrace data")
    parser.add_argument('-f', required=True,
        help="Base folder containing cleaned SELEX sequences")
    parser.add_argument('-i', required=True,
        help="Base folder containing SelexTrace run")
    parser.add_argument('-r', required=True, type=int,
        help="Round this input file comes from")
    parser.add_argument('--isc', type=int, default=80,
        help="Score cutoff for Infernal. (Default 80)")
    parser.add_argument('-c', type=int, default=1,
        help="Number of CPUs to use (Default 1)")

    args = parser.parse_args()
    if args.c < 1:
        raise ValueError("ERROR: CPU count must be at least 1!")
    if args.isc < 0:
        raise ValueError("ERROR: RNAforester score cutoff must be >0!")

    basefolder = args.f
    if basefolder[-1] != "/":
        basefolder += "/"
    outfolder = args.i
    if outfolder[-1] != "/":
        outfolder += "/"

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
        if group == "fasta_groups" or group == "fasta_families":
            continue
        run_infernal("%s%s/cmfile.cm" % (outfolder, group), args.r, seqs,
                     outfolder + group + "/", cpus=args.c, score=args.isc)
    del seqs
    # get families and write out families list
    families = create_families(outfolder, args.r, outfile="family_matrix.txt")
    print len(groups), "groups, ", len(families), "families"
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

    print "Runtime: %0.2f hrs" % ((time() - secs) / 3600)