from sys import argv
from os import walk
from os.path import exists
from os import walk
from numpy import empty

from cogent.parse.fasta import MinimalFastaParser

#overlap_check.py basefolder #round

if __name__ == "__main__":
    basefolder = argv[1]
    if basefolder[:-1] != "/":
        basefolder += "/"

    rnd = "R" + argv[2]


    #create csv file for data
    overlapfile = open(basefolder + rnd + "overlap.csv", 'w')
    #empty first cell of header
    overlapfile.write(",")

    #get all group names
    groups = walk(basefolder).next()[1]
    if "fasta_groups" in groups:
        groups.remove("fasta_groups")
    sizegroups = len(groups)
    #get total number fo groups in the round, add 1 for looping
    print str(sizegroups) + " groups"
    #write out all group names for column
    for x in range(sizegroups):
        overlapfile.write(groups[x] + ",")
    overlapfile.write("\n")
    #compare all group sequences to other sequences
    countmatrix = empty(shape=(sizegroups, sizegroups), dtype=int)
    countmatrix.fill(-1)
    for group, gfolder in enumerate(groups):
        #write out groupname to csv file as row indicator
        overlapfile.write(gfolder + ",")

        #compare seq counts over all groups
        for x in range(group):
            overlapfile.write("-,")
        #load in headers of current group
        headers = set([])
        firstgroup = open(basefolder + gfolder + "/" + rnd + "hits.fna")
        #get rid of header
        for header, seq in MinimalFastaParser(firstgroup):
            headers.add(header.split(" ")[0])
        firstgroup.close()
        dupefile = open(basefolder + gfolder + "/" + rnd + "dupes.txt", 'w')
        #go through all other groups and compare sequence headers
        for secgroup in range(group, sizegroups):
            if not exists(basefolder + groups[secgroup] + "/" + rnd + "hits.fna"):
                countmatrix[group][secgroup] = -1
                overlapfile.write("-,")
                continue
            count = 0
            if gfolder != groups[secgroup]:  # only do comparison if needed
                compfile = open(basefolder + groups[secgroup] + "/" + rnd + "hits.fna")
                dupefile.write(groups[secgroup] + "\n")
                for header, seq in MinimalFastaParser(compfile):
                    if header.split(" ")[0] in headers:
                        dupefile.write("%s\n%s\n" % (header, seq))
                        count += 1
                compfile.close()
            else:  # we are comparing group to itself so only need seq count
                count = len(headers)
            countmatrix[group][secgroup] = count
            overlapfile.write(str(count) + ",")
        overlapfile.write("\n")
    overlapfile.write("\n")

    #now create percentages matrix
    for group in range(sizegroups):
        overlapfile.write(groups[group] + ",")
        if countmatrix[group][group] != -1:
             #write out the dashes needed to line up the matrix in file
            for x in range(group):
                overlapfile.write("-,")

            groupcount = float(countmatrix[group][group])
            for secgroup in range(group, sizegroups):
                print "SECGROUP: ", countmatrix[group][secgroup]
                if countmatrix[group][secgroup] != -1:
                    print countmatrix[group][secgroup]/groupcount*100
                    overlapfile.write(str(countmatrix[group][secgroup]/groupcount*100) + ",")
                else:
                    overlapfile.write("-,")
        else:
            for secgroup in range(sizegroups):
                overlapfile.write("-,")
        overlapfile.write("\n")
