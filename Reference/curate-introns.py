#!/usr/bin/env python

#Jill E. Moore - Jill.Elizabeth.Moore@gmail.com
#Weng Lab - UMass Medical School
#Moore Lab Collaboration
#SPombe Splicing Project
#Updated June 2017

#python curate-introns.py Exons.bed > Introns.bed


import sys, numpy

def Create_Transcript_Dict(exonBed):
    tDict={}
    for line in exonBed:
        line=line.rstrip().split("\t")
        if line[3] not in tDict:
            tDict[line[3]]={"chrom":line[0],"gene":line[6],"strand":line[5], \
                            "exons":[[int(line[1]),int(line[2])]]}
        else:
            tDict[line[3]]["exons"].append([int(line[1]),int(line[2])])
    return tDict

def Print_Introns(transcript, transEntry):
    gene=transEntry["gene"]
    strand=transEntry["strand"]
    chrom=transEntry["chrom"]
    exons=numpy.array(transEntry["exons"])
    exons=exons[numpy.argsort(exons[:,0])]
    if strand == "+":
        start=exons[0][1]
        for i in range(1,len(exons)):
            print chrom+"\t"+str(start+1)+"\t"+str(exons[i][0]-1)+"\t"+ \
            transcript.rstrip()+"_intron_"+str(i)+"\t.\t"+strand+ \
            "\t"+ transcript.rstrip()+"_intron_"+str(i)+"\t" + \
            gene.rstrip()+"_intron_"+str(i)
            start=exons[i][1]
    elif strand == "-":
        start=exons[0][1]
        for i in range(1,len(exons)):
            print chrom+"\t"+str(start+1)+"\t"+str(exons[i][0]-1)+"\t"+\
            transcript.rstrip()+"_intron_"+str(len(exons)-i)+"\t.\t"+strand+ \
            "\t"+ transcript.rstrip()+"_intron_"+str(len(exons)-i)+"\t" + \
            gene.rstrip()+"_intron_"+str(len(exons)-i)
            start=exons[i][1]

exonBed=open(sys.argv[1])
tDict=Create_Transcript_Dict(exonBed)

for transcript in tDict:
    if len(tDict[transcript]["exons"]) > 1:
        Print_Introns(transcript, tDict[transcript])

