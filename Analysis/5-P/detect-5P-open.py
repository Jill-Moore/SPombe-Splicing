#Jill E Moore
#Weng Lab 
#UMass Medical School
#Updated Aug 2016

#Script for determining statistical significance of 5' phosphate pileups at position +1 of s pombe introns
#Run as:
#python detect-5P-5END.py regions.bed data.bam > data.output

import sys, pysam, re, numpy, scipy.stats
from scipy.stats import binom_test
from math import erf,sqrt


def Retrieve_Reads(bam, chrom, start, stop, strand):
    readCount=0
    for entry in bam.fetch(chrom, start-1, stop-1):
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
	    if strand == "+" and int(entry[3])+1 == start and not x.is_reverse:
                        readCount +=1
            elif strand == "-" and x.is_reverse and int(entry[3])+readLength == start:
                        readCount +=1
    return readCount

def Determine_Significance(bam, chrom, start, stop, strand):
    nullArray=[]
    if strand == "+":
	if start-12 > 1:
	    for i in range(start-12,start+13):
                nullArray.append(Retrieve_Reads(bam, chrom, i, i+1, strand))
	    fivePrime=Retrieve_Reads(bam, chrom, start, start+1, strand)
	else:
	    for i in range(1,start+25):
                nullArray.append(Retrieve_Reads(bam, chrom, i, i+1, strand))
	    fivePrime=Retrieve_Reads(bam, chrom, start, start+1, strand)
    elif strand == "-":
	if stop-12 > 1:
	    for i in range(stop-12,stop+13):
                nullArray.append(Retrieve_Reads(bam, chrom, i, i+1, strand))
	    fivePrime=Retrieve_Reads(bam, chrom, stop, stop+1, strand)        
	else:
	    for i in range(1,stop+25):
                nullArray.append(Retrieve_Reads(bam, chrom, i, i+1, strand))
	    fivePrime=Retrieve_Reads(bam, chrom, stop, stop+1, strand)
    totalReads=sum(nullArray)
    if totalReads == 0:
	p=1
    else:
	p=binom_test(fivePrime,totalReads,1/25.0)
    return fivePrime, totalReads, p
        
    
introns=open(sys.argv[1])
bam=pysam.Samfile(sys.argv[2], "rb")

for intron in introns:
    intron=intron.rstrip().split("\t")
    for i in range(int(intron[1]),int(intron[2])+1):
	readCount, totalReads, p = Determine_Significance(bam, intron[0], i, i, intron[5])
	if totalReads > 0:
    	    print intron[0]+"\t"+str(i)+"\t"+intron[3]+"\t"+intron[5]+"\t"+ str(readCount) +"\t" + str(totalReads)+"\t"+str(p)+ "\t" + str(readCount/float(totalReads))
	else:
	    print intron[0]+"\t"+str(i)+"\t"+intron[3]+"\t"+intron[5]+"\t"+ str(readCount) +"\t" + str(totalReads)+"\t"+str(p)+ "\t" + "0"

        
