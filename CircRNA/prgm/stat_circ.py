#!/usr/bin/env python
# coding=utf8

##Import modules
import os
import sys
import re
import urllib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
import subprocess

#import Bio.Graphics.GenomeDiagram

if __name__ == "__main__":
    file1=sys.argv[1]
    f1=open(file1,'r')
    b2=True
    f2=open(sys.argv[2],'w')
    gbk_ref=sys.argv[3]
    gbk_new=sys.argv[4]
    genome_size=int(sys.argv[5])
  
    dicmatch={} 
    dictotal={}
    dicstart={}
    dicend={}
    ### read the bam file and add each match in a dictionnary dictotal
    ### if dictotal contains already a match with this id, then the junction supported by the
    ### read is added in the dictionnary dicmatch    
    for line in f1:
        matchinit=re.search(r'^([^\s]*)\s+(\d+)\s+.*\s+(\d+)\s+\d+\s+.*\s+\d+\s+\d+\s+([A-Z]+).*',line)
	start=int(matchinit.group(3))
        length=len(matchinit.group(4))
        end=int(start)+length-1
        id =matchinit.group(1)
        
        if not dictotal.has_key(id):
            dictotal[id]=(start,end)
        else :
            if start>dictotal[id][0]:
                startm=dictotal[id][0]
                endm=end
            else :
                startm=start
                endm=dictotal[id][1]
            if not dicmatch.has_key(startm):
                dicmatch[startm]={}
            if not dicmatch[startm].has_key(endm):
                dicmatch[startm][endm]=0
            dicmatch[startm][endm]=1+dicmatch[startm][endm]

    f1.close()
    ### dicmatch contains the weighted list of end for each start
    ### Now for each junctions in dicmatch we search its annotations in the genbank file
    diclocus={}
    record=SeqIO.read(gbk_ref,"genbank")
    record_new=SeqIO.read(gbk_new,"genbank")
    new_locus=len(record.features)-len(SeqIO.read(gbk_ref,"genbank").features)+1
    new_feature=new_locus

    mem=genome_size
    for k in sorted(dicmatch):
        for k1, v1 in dicmatch[k].items():
	    feature_match=0
            feature_m=SeqFeature(FeatureLocation(0,genome_size), type="misc_RNA")
	    mem=abs(feature_m.location.nofuzzy_start-k)+abs(feature_m.location.nofuzzy_end-k1)
            for feature in record.features:
		mem_start=feature.location.nofuzzy_start
		mem_end=feature.location.nofuzzy_end
                if ( ((k1+k)/2 in feature) and (feature.location.nofuzzy_start>0 or feature.location.nofuzzy_end < genome_size) ) :
                    feature_match=1
                    if not "locus_tag" in feature.qualifiers :
                        feature.qualifiers["locus_tag"]=[]
                        feature.qualifiers["locus_tag"].append("NEW_Locus"+str(new_locus))
                        new_locus=new_locus+1
       		    if (abs(feature.location.nofuzzy_start-k)+abs(feature.location.nofuzzy_end-k1)<mem):
                        feature_m=feature
                        mem=abs(feature.location.nofuzzy_start-k)+abs(feature.location.nofuzzy_end-k1) 
			
            if (feature_match==0 or feature_m.location.nofuzzy_end-feature_m.location.nofuzzy_start>3*(k1-k) ):
                feature_new=SeqFeature(FeatureLocation(k1,k), type="misc_RNA")
		feature_new.qualifiers["locus_tag"]=[]
                if (feature_match>0):
                    feature_new.qualifiers["locus_tag"].append(feature_m.qualifiers["locus_tag"][0]+"_NA"+str(new_feature))
                else :
                    feature_new.qualifiers["locus_tag"].append("NA"+str(new_feature))
                new_feature=new_feature+1
                record_new.features.append(feature_new)
		feature_m=feature_new
            
	    if not diclocus.has_key(feature_m.qualifiers["locus_tag"][0]):
                diclocus[feature_m.qualifiers["locus_tag"][0]]={}
            if not diclocus[feature_m.qualifiers["locus_tag"][0]].has_key(k):
                diclocus[feature_m.qualifiers["locus_tag"][0]][k]=[]
            if diclocus[feature_m.qualifiers["locus_tag"][0]][k].count(k1)==0:
                diclocus[feature_m.qualifiers["locus_tag"][0]][k].append(k1)
            
    print("genbank finished\n")		   
    SeqIO.write(record_new,gbk_new,"genbank")


    for locus, ensstart in diclocus.items():
        nb_junc=0
        nb_read=0
        junc_mem=0
        junc=0
        start_junc=0
        end_junc=0
        for start, endlist in ensstart.items() :
            for end in endlist :
                nb_junc=nb_junc+1
                j=3
                junc=dicmatch[start][end]
                while (j>0) :
                    j=j-1
                    i=3
                    while (i+j>1 and i>0):
                        i=i-1
                        if diclocus[locus].has_key(start-j):
                            if end-i in diclocus[locus][start-j]:
                                junc=junc+dicmatch[start-j][end-i]
                            if end+i in diclocus[locus][start-j] and i>0:
                                junc=junc+dicmatch[start-j][end+i]
                        if diclocus[locus].has_key(start+j) and j>0:
                            if end-i in diclocus[locus][start+j] :
                                junc=junc+dicmatch[start+j][end-i]
                            if end+i in diclocus[locus][start+j] and i>0:
                                junc=junc+dicmatch[start+j][end+i]
                if (junc>junc_mem):
                    junc_mem=junc
                    start_junc=start
                    end_junc=end

        	string=locus+" "+str(start)+", "+str(end)+" "+str(dicmatch[start][end])+" "+str(junc)+" "+str(nb_junc)+"\n"
        	f2.write(string)

    print "nombre de reads circulaires "+str(len(dictotal))
    print "nombre de locus ayant au moins un read circulaire"+str(len(diclocus))
