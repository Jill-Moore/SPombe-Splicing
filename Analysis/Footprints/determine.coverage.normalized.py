#Jill Moore
#UMass Medical School
#Weng Lab
#Updated June 2017

#Script to generate coverage of spliceosome footprint reads around 5' splice
#sites, 3' splice sites and branch points. Read counts are normalized for
#by intron using the total number of reads in the window

#Usage: python determine.coverage.normalized.py bed bam

import sys, pysam, re, numpy

def Retrieve_Reads(bam, chrom, start, stop, strand, s1, s2):
    split=0
    exon=0
    intron=0
    across=0
    for entry in bam.fetch(chrom, start-1, stop-1):
        if (strand == "+" and not entry.is_reverse) or \
        (strand == "-" and entry.is_reverse):
            x=entry
            entry=str(entry).split()
            if "I" not in entry[5]:
                N= re.findall(r'\d+', entry[5])
                readLength= sum([ int(a) for a in N ])
            elif "I" in entry[5]:
                I=0
                S=""
                p=entry[5].split("I")
                for i in range(0,len(p)-1):
                    I += int(p[i][-1])
                    S += p[i][:-1]
                N= re.findall(r'\d+', entry[5])
                readLength= sum([ int(a) for a in N ])-I
            if "N" in entry[5]:
                N= re.findall(r'\d+', entry[5])
                if len(N) == 3:
                    if start <= int(entry[3])+int(N[0]) or \
                    start >= int(entry[3])+int(N[0])+int(N[1])+1:
                        split += 1
                else:
                    pass
            elif int(entry[3])+1 >= s1 and int(entry[3])+readLength <= s2:
                intron +=1
            elif  int(entry[3])+readLength <= s1-1 or int(entry[3])+1 >= s2+1:
                exon += 1
            else:
                across += 1
    return [split, exon, intron, across]

masterList=[]
for i in range(0,176):
    masterList.append([[],[],[],[]])

introns=open(sys.argv[1])
bam=pysam.Samfile(sys.argv[2], "rb")

for intron in introns:
    intron=intron.rstrip().split("\t")
    intronArray=[]
    t=0
    for i in range(int(intron[1]),int(intron[2])+1):
        x=Retrieve_Reads(bam, intron[0], i, i+1, intron[5], \
                                          int(intron[6]), int(intron[7]))
	t+=sum(x)
	intronArray.append(x)
    total=float(t)
    if intron[5] == "-":
        intronArray=intronArray[::-1]
    k=0
    for entry in intronArray:
	if total > 0:
            masterList[k][0].append(entry[0]/total)
            masterList[k][1].append(entry[1]/total)
            masterList[k][2].append(entry[2]/total)
            masterList[k][3].append(entry[3]/total)
	else:
	    masterList[k][0].append(0)
	    masterList[k][1].append(0)
	    masterList[k][2].append(0)
	    masterList[k][3].append(0)
        k+=1
for entry in masterList:
    print numpy.mean(entry[0]), "\t", numpy.mean(entry[1]), "\t", numpy.mean(entry[2]), "\t", numpy.mean(entry[3])

