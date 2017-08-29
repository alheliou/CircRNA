#!/usr/bin/env python
# coding=utf8

# attention dans cette version, utilisée depuis le 10 novembre, les reads circulaires sont comptés qu'ils n'aient que un ou deux matchs dans la région considérée.

##Import modules
import os
import sys
import re
import urllib
from Bio import SeqIO
#import matplotlib.pyplot as plt
import numpy
import subprocess
#import Bio.Graphics.GenomeDiagram

def analyseHB(outputpref,pos, nb_read_lin, nb_read_circ_entier, nb_read_circ,read_3,start, end, junction_id,seq_id,string):

    commande="samtools view -c "+outputpref+"_outl.sorted.bam "+seq_id+":"+str(start)+"-"+str(end)
    proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    (out, err) = proc.communicate()
    lin=int(out[:-1])
    print commande
    print lin 
    commande="samtools view "+outputpref+"_outc.sorted.bam "+seq_id+":"+str(start)+"-"+str(end)+"> temp.txt"
    proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    (out, err) = proc.communicate()
    print commande
    proc=subprocess.Popen("wc -l temp.txt",shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    (out, err) = proc.communicate() 	   
    file=open("temp.txt",'r')
    diccirc={}
    dic_junc={}

    circ_entier=0
    read_3_exp=0
    for line in file:
	matchinit=re.search(r'^.*:(\d+:\d+)\s+(\d+)\s+.*\s+(\d+)\s+\d+\s+.*\s+\d+\s+\d+\s+([A-Z]*).*',line)
        id =matchinit.group(1)
        
        
        if not diccirc.has_key(id):
            diccirc[id]=0
            dic_junc[id]=int(matchinit.group(3))
        else :
            diccirc[id]=diccirc[id]+1
            if diccirc[id]==1 :
                circ_entier+=1
                if (abs(dic_junc[id]-start)<4 and abs(int(matchinit.group(3))+len(matchinit.group(4))-end)<4):
                    read_3_exp+=1
            if diccirc[id]>1:
                print "pb avec"+id

    circ=len(diccirc)
    tot=circ+lin
    count=0
    if (tot >0) :
        tot=float(tot) ##test 22fev
        string=string+str(pos)+"\t"+junction_id+"\_"+outputpref+"\t"+str(start)+"\t"+str(end)+"\t"+str(100*read_3_exp/tot)+"\t"+str(format(100*circ_entier/tot-100*read_3_exp/tot,'.5f'))+"\t"+str(format(100-100*lin/tot-100*read_3_exp/tot,'.5f'))+"\t"+str(100*lin/tot)+"\t"+str(tot)+"\n"
        count=1
    nb_read_lin+=lin
    nb_read_circ+=circ
    nb_read_circ_entier+=circ_entier
    file.close()
    

    return (count,nb_read_lin,nb_read_circ_entier, nb_read_circ,read_3+read_3_exp,string)

def analyseHBlocus(outputpref,pos, nb_read_lin, nb_read_circ_entier, nb_read_circ,read_3, start, end,junc_start,junc_end, junction_id,seq_id,string):
    
    commande="samtools view -c "+outputpref+"_outl.sorted "+seq_id+":"+str(start)+"-"+str(end)
    proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    (out, err) = proc.communicate()
    lin=int(out[:-1])
    
    commande="samtools view "+outputpref+"_outc.sorted.bam "+seq_id+":"+str(start)+"-"+str(end)+"> temp.txt"
    proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    (out, err) = proc.communicate()
    
    file=open("temp.txt",'r')
    diccirc={}
    dic_junc={}
    
    circ_entier=0
    read_3_exp=0
    for line in file:
        matchinit=re.search(r'^.*:(\d+:\d+)\s+(\d+)\s+.*\s+(\d+)\s+\d+\s+.*\s+\d+\s+\d+\s+([A-Z]*).*',line)
	id =matchinit.group(1)
        
        
        if not diccirc.has_key(id):
            diccirc[id]=0
            dic_junc[id]=int(matchinit.group(3))
        else :
            diccirc[id]=diccirc[id]+1
            if diccirc[id]==1 :
                circ_entier+=1
                if (abs(dic_junc[id]-junc_start)<4 and abs(int(matchinit.group(3))+len(matchinit.group(4))-junc_end)<4):
                    read_3_exp+=1
            if diccirc[id]>1:
                print "pb avec"+id

    circ=len(diccirc)
    tot=circ+lin
    count=0
    if (tot >0) :
        tot=float(tot) ##test 22fev
        string=string+str(pos)+"\t"+junction_id+"\_"+outputpref+"\t"+str(start)+"\t"+str(end)+"\t"+str(100*read_3_exp/tot)+"\t"+str(format(100*circ_entier/tot-100*read_3_exp/tot, '.5f'))+"\t"+str(format(100-100*lin/tot-100*read_3_exp/tot,'.5f'))+"\t"+str(100*lin/tot)+"\t"+str(tot)+"\n"
        count=1
    nb_read_lin+=lin
    nb_read_circ+=circ
    nb_read_circ_entier+=circ_entier
    file.close()
    
    
    return (count,nb_read_lin,nb_read_circ_entier, nb_read_circ,read_3+read_3_exp,string)

if __name__ == "__main__":
    file1=sys.argv[1]
    f1=open(file1,'r')
    fo=open(sys.argv[3],'w')
    outputpref=sys.argv[2]
    seq_id=sys.argv[4]
    ref=sys.argv[5]
    dicstart={}
    pos_sn=0
    pos_r=0
    pos_t=0
    pos_nf=0
    pos_pab=0
    pos_pablocus=0
    pos_tlocus=0
    for line in f1:
        match=re.search(r'^([a-zA-Z0-9]*.*[a-zA-Z0-9]*)\s([0-9]*),\s([0-9]*)\s([0-9]*)\s([0-9]*)\s([0-9]*)',line)
        if not match: continue
        locus=match.group(1)
        start=int(match.group(2))
        end=int(match.group(3))
        read_0=int(match.group(4))
        read_3=int(match.group(5))
        
        j=3
        B=True
        while (B and j>-3):
            j=j-1
            if B and dicstart.has_key(start-j):
                i=3
                while (B and i>-3):
                    i=i-1
                    if dicstart[start-j].has_key(end-i):
                        dicstart[start-j][end-i][locus]=[read_0,read_3]
                        B=False
        if B:
            if not dicstart.has_key(start):
                dicstart[start]={}
            if not dicstart[start].has_key(end):
                dicstart[start][end]={}
            dicstart[start][end][locus]=[read_0,read_3]

    nb_matches=0
    data_read=[]
    data_read3=[]
    data_percent=[]
    dicsn={}
    dicr={}
    dict={}
    dicNF={}
    dicPAB={}
    readsn=0
    readr=0
    readt=0
    readPAB=0
    readNF=0
    pos=0
    sn_circ=0
    sn_3=0
    sn_lin=0
    NF_circ=0
    NF_3=0
    NF_lin=0
    PAB_circ=0
    PAB_3=0
    PAB_lin=0
    r_circ=0
    r_3=0
    r_lin=0
    t_circ=0
    t_3=0
    t_lin=0
    seq=SeqIO.read(ref,"fasta").seq
    
    Set_start=set([57375, 65333, 67951, 235439, 258066, 401610, 675407,893397, 960309, 1042252, 1292355])
    for start in sorted(dicstart): # ensend in dicstart.items():
        for end, locuslist in dicstart[start].items() :
            nb_matches=nb_matches+1
            #string=str(start)+", "+str(end)+" : "
            B=True
            for locus, value in locuslist.items() :
                # string=string+locus+" "
                if (B and ((len(locus)>3 and locus[3]=='s') or len(set([start-2,start-1,start,start+1,start+2]) & Set_start)>0 ) ):
                    max=0
                    for locus1, value in locuslist.items() :
                        if (value[1]>max):
                            max=value[1]
                            locus_m=locus1
                    if not (len(locus)>3 and locus[3]=='s'):
                        if 57375 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR46"
                            locuslist[locus]=locuslist[locus_m]
                        if 65333 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR22"
                            locuslist[locus]=locuslist[locus_m]
                        if 67951 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR32"
                            locuslist[locus]=locuslist[locus_m]
                        if 235439 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR49"
                            locuslist[locus]=locuslist[locus_m]
                        if 258066 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR31"
                            locuslist[locus]=locuslist[locus_m]
                        if 401610 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR27"
                            locuslist[locus]=locuslist[locus_m]
                        if 675407 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR25"
                            locuslist[locus]=locuslist[locus_m]
                        if 893397 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR30"
                            locuslist[locus]=locuslist[locus_m]
                        if 960309 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR56"
                            locuslist[locus]=locuslist[locus_m]
                        if 1042252 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR53"
                            locuslist[locus]=locuslist[locus_m]
                        if 1292355 in set([start-2,start-1,start,start+1,start+2]):
                            locus="sR41"
                            locuslist[locus]=locuslist[locus_m]
                    dicsn[locus]=(start+end)/2
                    readsn=readsn+locuslist[locus_m][1]
                    read_0=locuslist[locus_m][0]
                    read_3=locuslist[locus_m][1]
                    locus_m=locus
                    B=False
                if B and (len(locus)>3 and locus[3]=='t'):
                    max=0
                    for locus1, value in locuslist.items() :
                        if (value[1]>max):
                            max=value[1]
                            locus_m=locus1
                    dict[locus]=(start+end)/2
                    readt=readt+locuslist[locus_m][1]
                    read_0=locuslist[locus_m][0]
                    read_3=locuslist[locus_m][1]
                    print start,end,read_0,read_3, locus_m, locus
		    locus_m=locus
                    B=False
                if B and (len(locus)>3 and locus[3]=='r'):
                    max=0
                    for locus, value in locuslist.items() :
                        if (value[1]>max):
                            max=value[1]
                            locus_m=locus
                    dicr[locus]=(start+end)/2
                    readr=readr+locuslist[locus_m][1]
                    read_0=locuslist[locus_m][0]
                    read_3=locuslist[locus_m][1]
                    locus_m=locus
                    B=False
            if B:
                for locus, value in locuslist.items() :
                    if locus[0]=='P':
                        max=0
                        for locus, value in locuslist.items() :
                            if (value[1]>max):
                                max=value[1]
                                locus_m=locus
                        dicPAB[locus]=(start+end)/2
                        readPAB=readPAB+locuslist[locus_m][1]
                        read_0=locuslist[locus_m][0]
                        read_3=locuslist[locus_m][1]
                        locus_m=locus
                        B=False
                        break
            if B :
                for locus, value in locuslist.items() :
                    max=0
                    for locus, value in locuslist.items() :
                        if (value[1]>max):
                            max=value[1]
                            locus_m=locus
                    dicNF[locus_m]=(start+end)/2
                    readNF=readNF+locuslist[locus_m][1]
                    read_0=locuslist[locus_m][0]
                    read_3=locuslist[locus_m][1]
                    break
            nb_read_lin=0
            nb_read_circ=0
            nb_read_circ_entier=0
            read_3_exp=0
            string=""
            pos1=0
            if dicsn.has_key(locus_m) :
                pos1=pos_sn
            if dict.has_key(locus_m) :
                pos1=pos_t
            if dicNF.has_key(locus_m) :
                pos1=pos_nf
            if dicPAB.has_key(locus_m) :
                pos1=pos_pab
            if dicr.has_key(locus_m) :
                pos1=pos_r
            count=0
            sum=0
            (count,nb_read_lin,nb_read_circ_entier,nb_read_circ,read_3_exp,string)=analyseHB(outputpref,pos1, nb_read_lin,nb_read_circ_entier, nb_read_circ,read_3_exp,start, end,locus_m.replace("_","\_"),seq_id,string)
            ##
	    print nb_read_circ, nb_read_lin #
            ##
            sum+=count

            if dicsn.has_key(locus_m) :
                pos_sn+=count+1
            if dict.has_key(locus_m) :
                pos_t+=1+count
            if dicNF.has_key(locus_m) :
                pos_nf+=1+count
            if dicPAB.has_key(locus_m) :
                pos_pab+=1+count
            if dicr.has_key(locus_m) :
                pos_r+=1+count
            
            if not dicr.has_key(locus_m):
                data_read.append(read_0)
                data_read3.append(read_3)
		if (read_3==0 or nb_read_circ==0):
			data_percent.append(0)
		else:
                	data_percent.append((100*read_3)/(nb_read_circ))

            if 100*read_3>=50*nb_read_circ_entier:
                    stringforfile=""+locus_m+" "+str(start)+" "+str(end)+" "+str(read_3)+" "+str(nb_read_circ_entier)+" "+str(nb_read_circ)+" "+str(nb_read_lin)+"\n"
		    fo.write(stringforfile)
		    if dicsn.has_key(locus_m):
                        print string
                        sn_circ=sn_circ+nb_read_circ
                        sn_3=sn_3+read_3
                        sn_lin=sn_lin+nb_read_lin
                        pos=pos+1
                        print locus_m, start, end, seq[start-20:start+1], seq[end: end+21]
                    if dicNF.has_key(locus_m):
                        NF_circ=NF_circ+nb_read_circ
                        NF_3=NF_3+read_3
                        NF_lin=NF_lin+nb_read_lin
                        print locus_m, start, end, seq[start-20:start+1], seq[end: end+21]
                    if dicPAB.has_key(locus_m):
                        PAB_circ=PAB_circ+nb_read_circ
                        PAB_3=PAB_3+read_3
                        PAB_lin=PAB_lin+nb_read_lin
                        print locus_m, start, end, seq[start-20:start+1], seq[end: end+21]
                    
                    if dict.has_key(locus_m):
                        t_circ=t_circ+nb_read_circ
                        t_3=t_3+read_3
                        t_lin=t_lin+nb_read_lin
                        print locus_m, start, end, seq[start-20:start+1], seq[end: end+21]
                                
                    if dicr.has_key(locus_m):
                        r_circ=r_circ+nb_read_circ
                        r_3=r_3+read_3
                        r_lin=r_lin+nb_read_lin

            else :
                if dicsn.has_key(locus_m):
                    if locus_m[3:10]=="snRNA32":
                        print read_3, nb_read_circ_entier, string
                    dicsn.pop(locus_m,None)
                    readsn=readsn-read_3
                    pos_sn-=1+count
                if dicr.has_key(locus_m):
                    dicr.pop(locus_m,None)
                    readr=readr-read_3
                    pos_r-=1+count
                if dict.has_key(locus_m):
                    dict.pop(locus_m,None)
                    readt=readt-read_3
                    pos_t-=1+count
                if dicPAB.has_key(locus_m):
                    dicPAB.pop(locus_m,None)
                    readPAB=readPAB-read_3
                    pos_pab-=1+count
                if dicNF.has_key(locus_m):
                    dicNF.pop(locus_m,None)
                    readNF=readNF-read_3
                    pos_nf-=1+count
                                                                                           
                                                                                           
    nb_read_lin=0
    nb_read_circ=0
    nb_read_circ_entier=0
    read_3_exp=0
    string=""
    pos=0
    print nb_matches

    print "nombre de locus sn "+str(len(dicsn))
    print "nombre de locus t "+str(len(dict))
    print "nombre de locus r "+str(len(dicr))
    print "nombre de locus new "+str(len(dicNF))
    print "nombre de locus autres "+str(len(dicPAB))

    nb_locus=len(dicsn)+len(dict)+len(dicr)+len(dicNF)+len(dicPAB)
    string="0/"+str(360*len(dicsn)/nb_locus)+"/"+str(100*len(dicsn)/nb_locus)+"/"+str(len(dicsn))+"/CDbox/50,"
    string=string+str(360*len(dicsn)/nb_locus)+"/"+str(360*(len(dicsn)+len(dicPAB))/nb_locus)+"/"+str(100*(len(dicPAB))/nb_locus)+"/"+str(len(dicPAB))+"/PAB/10,"
    string=string+str(360*(len(dicsn)+len(dicPAB))/nb_locus)+"/"+str(360*(len(dicsn)+len(dicPAB)+len(dicNF))/nb_locus)+"/"+str(100*(len(dicNF))/nb_locus)+"/"+str(len(dicNF))+"/non-annotated/20,"
    string=string+str(360*(len(dicsn)+len(dicPAB)+len(dicNF))/nb_locus)+"/"+str(360*(len(dicsn)+len(dicPAB)+len(dicNF)+len(dict))/nb_locus)+"/"+str(100*(len(dict))/nb_locus)+"/"+str(len(dict))+"/tRNA/30,"
    string=string+str(360*(len(dicsn)+len(dict)+len(dicPAB)+len(dicNF))/nb_locus)+"/360/"+str(100*(len(dicr))/nb_locus)+"/"+str(len(dicr))+"/rRNA/40"

    print "nombre de locus total "+str(nb_locus)

    print "nombre de reads soutenant à 3 près la jonction d'un locus sn "+str(readsn)
    print "nombre de reads soutenant à 3 près la jonction d'un locus t "+str(readt)
    print "nombre de reads soutenant à 3 près la jonction d'un locus r "+str(readr)
    print "nombre de reads soutenant à 3 près la jonction d'un locus new "+str(readNF)
    print "nombre de reads soutenant à 3 près la jonction d'un locus autre "+str(readPAB)

    nb_read=readsn+readt+readr+readNF+readPAB
    if (nb_read==0):
	nb_read=1
    string="0/"+str(360*readsn/nb_read)+"/"+str(100*readsn/nb_read)+"/"+str(readsn)+"/CDbox/50,"
    string=string+str(360*readsn/nb_read)+"/"+str(360*(readsn+readPAB)/nb_read)+"/"+str(100*(readPAB)/nb_read)+"/"+str(readPAB)+"/PAB/10,"
    string=string+str(360*(readsn+readPAB)/nb_read)+"/"+str(360*(readsn+readPAB+readNF)/nb_read)+"/"+str(100*(readNF)/nb_read)+"/"+str(readNF)+"/non-annotated/20,"
    string=string+str(360*(readsn+readPAB+readNF)/nb_read)+"/"+str(360*(readsn+readPAB+readNF+readt)/nb_read)+"/"+str(100*(readt)/nb_read)+"/"+str(readt)+"/tRNA/30,"
    string=string+str(360*(readsn+readPAB+readNF+readt)/nb_read)+"/360/"+str(100*(readr)/nb_read)+"/"+str(readr)+"/rRNA/40"

    for k,v in dicsn.items():
        print k

    f1.close()
    fo.close()
    print len(data_read)
    #plt.subplot(311)
    #plt.hist(data_read,bins=1000, log=True)
    print "nombre de reads soutenant une jonction"
    print "moyenne "+str(numpy.mean(data_read))+" variance "+str(numpy.std(data_read))+ " médianne "+str(numpy.median(data_read))+ " dernier vingtieme "+str(numpy.percentile(data_read,95))+ " premier vingtieme "+str(numpy.percentile(data_read3,5))
    #plt.subplot(312)
    #plt.hist(data_read3,bins=1000, log=False)
    print "nombre de reads soutenant une jonction à 2 nucleotides pres"
    print "moyenne "+str(numpy.mean(data_read3))+" variance "+str(numpy.std(data_read3))+ " médianne "+str(numpy.median(data_read3))+ " dernier vingtieme "+str(numpy.percentile(data_read3,95))+ " premier vingtieme "+str(numpy.percentile(data_read3,5))
    #plt.subplot(313)
    #plt.hist(data_percent,bins=1000,log=False)
    print " % des reads de ce locus soutenant une jonction"
    print "moyenne "+str(numpy.mean(data_percent))+" variance "+str(numpy.std(data_percent))+ " médianne "+str(numpy.median(data_percent))+ " dernier vingtieme "+str(numpy.percentile(data_percent,95))+ " premier vingtieme "+str(numpy.percentile(data_read3,5))


