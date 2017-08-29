#!/usr/bin/env python
# coding=utf8

##Import modules
import os
import sys
import re
import urllib
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline


def main(filexml, fileoutc, fileoutl,fileref, filefasta):
    reads=SeqIO.index(filefasta,"fasta") 
    ref=str(SeqIO.read(fileref,"fasta").seq)
    id=str(SeqIO.read(fileref,"fasta").id)
    print "reference: ",id, len(ref)
    file_in=open(filexml,"r")
    file_out=open(fileoutc,"w")
    file_outlin=open(fileoutl,"w")
    file_out.write("@HD\tVN:1.0\tS0:unsorted\n@SQ\tSN:"+id+"\tLN:1765118\n@PG\tID:blast_analysis.py\n")
    file_outlin.write("@HD\tVN:1.0\tS0:unsorted\n@SQ\tSN:"+id+"\tLN:1765118\n@PG\tID:blast_analysis.py\n")
    blast_records = NCBIXML.parse(file_in)
    
    
    for blast_record in blast_records:
        longueur_query=blast_record.query_letters
        for alignment in blast_record.alignments:
            b=0
            for hsp in alignment.hsps:
                if b>0:
                    break
                if hsp.query_start<hsp.query_end:
                    query_start=hsp.query_start
                    query_end=hsp.query_end
                    hspq=1
                else:
                    query_start=hsp.query_end
                    query_end=hsp.query_start
                    hspq=-1
                if hsp.sbjct_start<hsp.sbjct_end:
                    sbjct_start=hsp.sbjct_start
                    sbjct_end=hsp.sbjct_end
                    hsps=1
                else:
                    sbjct_start=hsp.sbjct_end
                    sbjct_end=hsp.sbjct_start
                    hsps=-1
                
                sens=hspq*hsps
                if len(hsp.query)==longueur_query:
                    b=5
                    if sens==1 :
                        sens="0"
                    else :
                        sens="16"
                    
                    (seq,cigar)=align(hsp.query,hsp.sbjct)
                    file_outlin.write(blast_record.query+"\t"+sens+"\t"+id+"\t"+str(sbjct_start)+"\t255\t"+cigar+"\t*\t0\t0\t"+seq+"\t*\n")
                    break
                    
                if (blast_record.query.split(" ")[0]=="read:3:8"):
		    print query_start, query_end, sbjct_start, sbjct_end 
                
                for hsp2 in alignment.hsps:
                    if hsp2.query_start<hsp2.query_end:
                        query2_start=hsp2.query_start
                        query2_end=hsp2.query_end
                        hsp2q=1
                    else:
                        query2_start=hsp2.query_end
                        query2_end=hsp2.query_start
                        hsp2q=-1
                    if hsp2.sbjct_start<hsp2.sbjct_end:
                        sbjct2_start=hsp2.sbjct_start
                        sbjct2_end=hsp2.sbjct_end
                        hsp2s=1
                    else:
                        sbjct2_start=hsp2.sbjct_end
                        sbjct2_end=hsp2.sbjct_start
                        hsp2s=-1

                    diff_start_query=query_start-query2_start
                    diff_end_query=query_end-query2_end
                    split_query=query_start-query2_end
                    split_sbjct=sbjct2_start-sbjct_end
                    diff_sbjct=sbjct2_start-sbjct_start
                    longueur=hsp.align_length+hsp2.align_length
                
                    if (query2_start<3 and query_end>longueur_query-2 and diff_start_query>0 and diff_end_query>0  and diff_sbjct*sens>0 and sens*split_sbjct>=0 and split_query<=2 and hspq==hsp2q and hsp2s==hsps and diff_sbjct<10000 and diff_sbjct>-10000 and b==0):
                        b=1
                        if sens==1 :
                            sens="0"
                        else :
                            sens="16"
                        
                        if (hsps==-1):
                            str_query=str(Seq(hsp.query,IUPAC.unambiguous_dna).reverse_complement())
                        else:
                            str_query=hsp.query
                        if (hsp2s==-1):
                            str2_query=str(Seq(hsp2.query,IUPAC.unambiguous_dna).reverse_complement())
                        else:
                            str2_query=hsp2.query
                        

                        (seq,cigar)=align(hsp.query,hsp.sbjct)
                        if (hsps==-1):
                            str_seq=str(Seq(seq,IUPAC.unambiguous_dna).reverse_complement())
                        else:
                            str_seq=seq
                        file_out.write(blast_record.query+"\t"+sens+"\t"+id+"\t"+str(sbjct_start)+"\t255\t"+cigar+"\t*\t0\t0\t"+str_seq+"\t*\n")
                        (seq,cigar)=align(hsp2.query,hsp2.sbjct)
                        if (hsp2s==-1):
                            str_seq=str(Seq(seq,IUPAC.unambiguous_dna).reverse_complement())
                        else:
                            str_seq=seq
                        file_out.write(blast_record.query+"\t"+sens+"\t"+id+"\t"+str(sbjct2_start)+"\t255\t"+cigar+"\t*\t0\t0\t"+str_seq+"\t*\n")
### partie à enlever quand on ne veut que les big matches
		if( 0==0):            
                    if (query_start==1 and longueur_query-query_end<11 and longueur_query-query_end>5 and b==0):
                        b=2
                        if hsps==1 :
                            if hspq==1 :
                            
                                b=motif_search(ref[sbjct_start-101:sbjct_start-1].upper(), str(reads[blast_record.query.split(" ")[0]].seq[query_end:longueur_query]).upper(),blast_record.query.split(" ")[0],1,id,sbjct_start-101,1,-1,file_out)
                        
                            if hspq==-1 :
                                b=motif_search(ref[sbjct_end:sbjct_end+100].upper(), str(reads[blast_record.query.split(" ")[0]].seq[query_end:longueur_query].reverse_complement()).upper(),blast_record.query.split(" ")[0],-1,id,sbjct_end,1,1,file_out)
                        if hsps==-1 :
                            if hspq==1 :
                                b=motif_search(str(Seq(ref[sbjct_end:sbjct_end+100],IUPAC.unambiguous_dna)).upper(), str(reads[blast_record.query.split(" ")[0]].seq[query_end:longueur_query].reverse_complement()).upper(),blast_record.query.split(" ")[0],-1,id,sbjct_end,1,1, file_out)
                            if hspq==-1 :
                                b=motif_search(str(Seq(ref[sbjct_start-101:sbjct_start-1],IUPAC.unambiguous_dna)).upper(), str(reads[blast_record.query.split(" ")[0]].seq[query_end:longueur_query]).upper(),blast_record.query.split(" ")[0],1,id,sbjct_start-101,1,-1,file_out)

                
                    if( query_end==longueur_query and query_start<=11 and query_start>6 and b==0):
                        b=2
                        if hsps==1 :
                            if hspq==1 :
                                b=motif_search(ref[sbjct_end:sbjct_end+100].upper(), str(reads[blast_record.query.split(" ")[0]].seq[0:query_start-1]).upper(),blast_record.query.split(" ")[0], 1, id, sbjct_end,1,1, file_out)
                            if hspq==-1 :
                                b=motif_search(ref[sbjct_start-101:sbjct_start-1].upper(), str(reads[blast_record.query.split(" ")[0]].seq[0:query_start-1].reverse_complement()).upper(),blast_record.query.split(" ")[0], -1, id, sbjct_start-101, 1, -1,file_out)
                        if hsps==-1 :
                            if hspq==1 :
                                b=motif_search(str(Seq(ref[sbjct_start-101:sbjct_start-1],IUPAC.unambiguous_dna)).upper(), str(reads[blast_record.query.split(" ")[0]].seq[0:query_start-1].reverse_complement()).upper(),blast_record.query.split(" ")[0], -1, id,sbjct_start-101,1,-1,file_out )
                            if hspq==-1 :
                                b=motif_search(str(Seq(ref[sbjct_end:sbjct_end+100],IUPAC.unambiguous_dna)).upper(), str(reads[blast_record.query.split(" ")[0]].seq[0:query_start-1]).upper(),blast_record.query.split(" ")[0], 1, id, sbjct_end, 1,1,file_out)
###
                if (b==3) :
                    if sens==1 :
                        sens="0"
                    else :
                        sens="16"
                    if (hsps==-1):
                        str_query=str(Seq(hsp.query,IUPAC.unambiguous_dna).reverse_complement())
                    else:
                        str_query=hsp.query
                    (seq,cigar)=align(hsp.query,hsp.sbjct)
                    if (hsps==-1):
                        str_seq=str(Seq(seq,IUPAC.unambiguous_dna).reverse_complement())
                    else:
                        str_seq=seq
                    file_out.write(blast_record.query.split(" ")[0]+"\t"+sens+"\t"+id+"\t"+str(sbjct_start)+"\t255\t"+cigar+"\t*\t0\t0\t"+str_seq+"\t*\n")
        
                if ( (b==0 or b==2)and query_end >longueur_query-3 and query_start<=3):
		    if sens==1 :
                        sens="0"
                    else :
                        sens="16"
                    if (hsps==-1):
                        str_query=str(Seq(hsp.query,IUPAC.unambiguous_dna).reverse_complement())
                    else:
                        str_query=hsp.query
                    (seq,cigar)=align(hsp.query,hsp.sbjct)
                    if (hsps==-1):
                        str_seq=str(Seq(seq,IUPAC.unambiguous_dna).reverse_complement())
                    else:
                        str_seq=seq
                    file_outlin.write(blast_record.query.split(" ")[0]+"\t"+sens+"\t"+id+"\t"+str(sbjct_start)+"\t255\t"+cigar+"\t*\t0\t0\t"+str_seq+"\t*\n")

    file_out.close()
    file_outlin.close()
    file_in.close()

def motif_search(ref, query,blast_record_query, sens, id, sbjct_pos, sbjct_sens,search_sens, file_out):
    b=0
    k=2
    
    if len(query)<8:
        k=1
    if len(query)<4:
        k=0
    if len(query)<3:
        return 2
    if search_sens==1:
        index=0
        while (index<=len(ref)-len(query)):
            score=0
            i=0
            cigar=""
            M=0
            S=0
            while (i<len(query) and score<=k):
                if ref[index+i]==query[i]:
                    if S!=0:
                        cigar=cigar+str(S)+"S"
                        S=0
                    M=M+1
                else :
                    if M!=0:
                        cigar=cigar+str(M)+"M"
                        M=0
                    score=score+1
                    S=S+1
                i=i+1
                if score>k:
                    break
            if (score<=k and i==len(query)):
                b=1
                break
            index=index+1
    
        sbjct_start=sbjct_pos+index+1

    if search_sens==-1:
        index=len(ref)-len(query)
        while (index>0):
            score=0
            i=0
            cigar=""
            M=0
            S=0
            while (i<len(query) and score<=k):
                if ref[index+i]==query[i]:
                    if S!=0:
                        cigar=cigar+str(S)+"S"
                        S=0
                    M=M+1
                else :
                    if M!=0:
                        cigar=cigar+str(M)+"M"
                        M=0
                    score=score+1
                    S=S+1
                i=i+1
                if score>k:
                    break
            if (score<=k and i==len(query)):
                b=1
                break
            index=index-1

        sbjct_start=sbjct_pos+index+1

    if (b==1) :
        if sens==1 :
            sens="0"
        else :
            sens="16"
        if M!=0:
            cigar=cigar+str(M)+"M"
        if S!=0:
            cigar=cigar+str(S)+"S"
        
        str_query=query
        
        
        file_out.write(blast_record_query+"\t"+sens+"\t"+id+"\t"+str(sbjct_start)+"\t255\t"+cigar+"\t*\t0\t0\t"+str_query+"\t*\n")
        return 3
            
    return 2

def align(query,sbjct):
    cigar=""
    seq=""
    if len(query)==len(sbjct) :
        index=0
        m=0
        i=0
        s=0
        d=0
        while index < len(query):
            if query[index]==sbjct[index]:
                m=m+1
                seq=seq+query[index]
                if i!=0:
                    cigar=cigar+str(i)+"I"
                    i=0
                if s!=0:
                    cigar=cigar+str(s)+"S"
                    s=0
                if d!=0:
                    cigar=cigar+str(d)+"D"
                    d=0
            if query[index]=="-":
                d=d+1
                if i!=0:
                    cigar=cigar+str(i)+"I"
                    i=0
                if s!=0:
                    cigar=cigar+str(s)+"S"
                    s=0
                if m!=0:
                    cigar=cigar+str(m)+"M"
                    m=0
            if sbjct[index]=="-":
                i=i+1
                if m!=0:
                    cigar=cigar+str(m)+"M"
                    m=0
                if s!=0:
                    cigar=cigar+str(s)+"S"
                    s=0
                if d!=0:
                    cigar=cigar+str(d)+"D"
                    d=0
                seq=seq+query[index]
            if (query[index]!="-" and sbjct[index]!="-" and query[index]!=sbjct[index]) :
                s=s+1
                if i!=0:
                    cigar=cigar+str(i)+"I"
                    i=0
                if m!=0:
                    cigar=cigar+str(m)+"M"
                    m=0
                if d!=0:
                    cigar=cigar+str(d)+"D"
                    d=0
                seq=seq+query[index]
            index=index+1
    else :
        print " not the same length\n"


    if i!=0:
        cigar=cigar+str(i)+"I"
        i=0
    if s!=0:
        cigar=cigar+str(s)+"S"
        s=0
    if d!=0:
        cigar=cigar+str(d)+"D"
        d=0
    if m!=0:
        cigar=cigar+str(m)+"M"
        m=0

    return (seq,cigar)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
