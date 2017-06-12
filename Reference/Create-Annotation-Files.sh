#!/bin/bash 

#Jill E. Moore - Jill.Elizabeth.Moore@gmail.com
#Weng Lab - UMass Medical School
#Moore Lab Collaboration
#SPombe Splicing Project
#Updated June 2017

#Script which takes input gtf file and creates annotation files
#We used Schizosaccharomyces_pombe.ASM294v2.35.gtf downloaded 06-08-2017
#from ensembl for our analyses

#./Create-Annotation-Files.sh SPombe.gtf

gtf=$1

awk '{if ($3== "gene") print $1 "\t" $4 "\t" $5 "\t" $10 "\t" "." "\t" $7 "\t" \
     $12}' $gtf | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' \
      > Genes.bed

awk '{if ($3== "transcript") print $1 "\t" $4 "\t" $5 "\t" $12 "\t" "." "\t" $7 \
     "\t" $10}' $gtf | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' \
     > Transcripts.bed

awk '{if ($3== "exon") print $1 "\t" $4 "\t" $5 "\t" $12 "\t" "." "\t" $7 "\t" \
     $10}' $gtf | awk '{gsub(/;/," ");print}' | awk '{gsub(/"/,"");print}' \
     > Exons.bed

python curate-introns.py Exons.bed > Introns.bed

awk '{print $1 "\t" $2 "\t" $3}' Genes.bed | sort -k1,1 -k2,2n > sorted
bedtools merge -i sorted > merged

python curate-intergenic-regions.py merged
rm merged sorted
