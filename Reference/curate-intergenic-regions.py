#!/usr/bin/env python

#Jill E. Moore - Jill.Elizabeth.Moore@gmail.com
#Weng Lab - UMass Medical School
#Moore Lab Collaboration
#SPombe Splicing Project
#Updated June 2017

#python curate-introns.py Exons.bed > Introns.bed

import sys, numpy

def Create_Region_Dict(genes):
    regionDict={}
    for line in genes:
        line=line.rstrip().split("\t")
        if line[0] not in regionDict:
            regionDict[line[0]]=[[int(line[1]),int(line[2])]]
        else:
            regionDict[line[0]].append([int(line[1]),int(line[2])])
    return regionDict

def Print_Intergenic_Regions(regionDict):
    watsonOut=open("Intergenic-Regions-Watson.bed", "w+")
    crickOut=open("Intergenic-Regions-Crick.bed", "w+")
    for chrom in regionDict:
        regions=numpy.array(regionDict[chrom])
        regions=regions[numpy.argsort(regions[:,0])]
        start=regions[0][1]
        for i in range(1,len(regions)):
            print >> watsonOut,chrom+"\t"+str(start+1)+"\t"+str(regions[i][0]-1) \
            +"\t"+ "Intergenic_Region_"+str(i)+"\t.\t+"
            print >> crickOut,chrom+"\t"+str(start+1)+"\t"+str(regions[i][0]-1) \
            +"\t"+ "Intergenic_Region_"+str(i)+"\t.\t-"
            start=regions[i][1]
    watsonOut.close()
    crickOut.close()

genes=open(sys.argv[1])
regionDict=Create_Region_Dict(genes)
Print_Intergenic_Regions(regionDict)




