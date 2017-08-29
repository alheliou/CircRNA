#!/usr/bin/env python
# coding=utf8

# afficher la distribution en nombre et en pourcentage, du nombre de reads qui soutiennent une jonction
# comparer aussi les positions de jonctions en fonction de SRR et HB
# refaire la prediction de structure avec RNAfold,(mfe) pour les sno et le 5S
# est-ce que les sno ont une fonctions chez Paby -> est-ce que l'on peut trouver l'ARN complémentaire (5 nucleotides avant la boite D), de taille 6-8 ou plus


# faire les diagrames avec les manip traitées RNAse, regarder qu'elles sont les jonctions qui restent

# C/D vbox sRNA-guided 2'-O-methylation patterns of archaela rRNA molecules
# bacterie : aquifex aeolicus


## attention nouveauté du 25 novembre, pour récupérer le snRNA32 avec une taille correcte, ajout de end_junc-start_junc<end-start qui en cas d'egalité parfaite entre deux jonctions, donne la préférence à celle qui est plus grande. Cela n'influence véritablement que le snRNA32, pour les autres locus, la différence ne se joue qu'a quelques nucléotides


##Import modules
import os
import sys
import re
import urllib
from Bio import SeqIO
#import matplotlib.pyplot as plt
import numpy
#import Bio.Graphics.GenomeDiagram

if __name__ == "__main__":
    file1=sys.argv[1]
    f1=open(file1,'r')
    files=[f1]
    read_mini=1
    file_mini=0
    if len(sys.argv)> 2:
	read_mini=3
	file_mini=1
        f2=open(sys.argv[2],'r')
        files.append(f2)
        if len(sys.argv) > 3:
            f3=open(sys.argv[3],'r')
            files.append(f3)
            if len(sys.argv)>4:
                f4=open(sys.argv[4],'r')
                files.append(f4)
                if len(sys.argv)>5:
                    f5=open(sys.argv[5],'r')
                    files.append(f5)
    dict={}
    dicsn={}
    dicr={}
    dicnew={}
    dicPAB={}
    diclocus={}
    dicmatch={}
    exp_int=-1
    #R=[164349,265669,35811,234948,37113394]
    #C=[36692,41952,12784,60302,4293888]


### lecture du fichier, a chauqe locus on associe les jonctions avec les infos : read_0, read_3 et nb de file
    for f in files:
        for line in f:
            match=re.search(r'^([a-zA-Z0-9]*.*[a-zA-Z0-9]*)\s([0-9]*),\s([0-9]*)\s([0-9]*)\s([0-9]*)\s([0-9]*)',line)
            if not match :
                print "pb with "+line
                continue
            locus=match.group(1)
            start=int(match.group(2))
            end=int(match.group(3))
            read_0=int(match.group(4))
            read_3=int(match.group(5))
            nb_junc=int(match.group(6))
            if not diclocus.has_key(locus):
                diclocus[locus]={}
            if not diclocus[locus].has_key(start):
                diclocus[locus][start]={}
            if not diclocus[locus][start].has_key(end):
                diclocus[locus][start][end]=[read_0,read_3,1]
            else :
                diclocus[locus][start][end][0]=diclocus[locus][start][end][0]+read_0
                diclocus[locus][start][end][1]=diclocus[locus][start][end][1]+read_3
                diclocus[locus][start][end][2]=diclocus[locus][start][end][2]+1
    f1.close()
    if len(sys.argv)> 2:
        f2.close()
    if len(sys.argv)> 3:
        f3.close()
    if len(sys.argv)> 4:
        f4.close()
    if len(sys.argv)> 5:
        f5.close()

## pour chaque locus, nous cherchons une jonction "identifiable"
    for locus, ensstart in diclocus.items():
        nb_junc=0
        junc_mem=0
        junc0_mem=0
        junc=0
        start_junc=0
        end_junc=0
        nb_file_max=1

        for start, endlist in ensstart.items() :
            for end in endlist :
                read_0= diclocus[locus][start][end][0]
                nb_file=diclocus[locus][start][end][2]
                j=3
                
                read_3=read_0
                while (j>0) :
                    j=j-1
                    i=3
                    while (i>0 and i+j>1):
                        i=i-1
                        if diclocus[locus].has_key(start-j):
                            if diclocus[locus][start-j].has_key(end-i):
                                read_3=read_3+diclocus[locus][start-j][end-i][0]
                            if diclocus[locus][start-j].has_key(end+i)and i>0:
                                read_3=read_3+diclocus[locus][start-j][end+i][0]
                        if diclocus[locus].has_key(start+j) and j>0:
                            if diclocus[locus][start+j].has_key(end-i):
                                read_3=read_3+diclocus[locus][start+j][end-i][0]
                            if diclocus[locus][start+j].has_key(end+i) and i>0:
                                read_3=read_3+diclocus[locus][start+j][end+i][0]
                if read_3>diclocus[locus][start][end][1] and diclocus[locus][start][end][2]==1 :
                    print read_3, diclocus[locus][start][end][1]
                    nb_file=2
                else :
                    nb_file=diclocus[locus][start][end][2]
                if ((read_3>junc_mem or (read_3==junc_mem and (read_0>junc0_mem or (read_0==junc0_mem and end_junc-start_junc<end-start)) ) ) and (nb_file==nb_file_max or nb_file >1)):
                    start_junc=start
                    end_junc=end
                    junc_mem=read_3
                    junc0_mem=read_0
                    nb_file_junc=nb_file
                    nb_file_max=nb_file

        read_0= diclocus[locus][start_junc][end_junc][0]
        read_3= junc_mem
        if read_3>0 :
            nb_file= nb_file_junc
            if (read_3>read_mini and nb_file>file_mini ):
                if len(locus)>3 and locus[3]=='s':
                    dicsn[locus]=1
                    print locus+" "+str(start_junc)+", "+str(end_junc)+" "+str(read_0)+" "+str(read_3)+" "+str(nb_junc)+" "+str(diclocus[locus][start_junc][end_junc][2])+" "+str(nb_file)
                elif len(locus)>3 and locus[3]=='t':
                    dict[locus]=1
                    print locus+" "+str(start_junc)+", "+str(end_junc)+" "+str(read_0)+" "+str(read_3)+" "+str(nb_junc)+" "+str(diclocus[locus][start_junc][end_junc][2])+" "+str(nb_file)
                elif len(locus)>3 and locus[3]=='r':
                    dicr[locus]=1
                    print locus+" "+str(start_junc)+", "+str(end_junc)+" "+str(read_0)+" "+str(read_3)+" "+str(nb_junc)+" "+str(diclocus[locus][start_junc][end_junc][2])+" "+str(nb_file)
                elif locus[0]=='N':
                    dicnew[locus]=1
                    print locus+" "+str(start_junc)+", "+str(end_junc)+" "+str(read_0)+" "+str(read_3)+" "+str(nb_junc)+" "+str(diclocus[locus][start_junc][end_junc][2])+" "+str(nb_file)
                else :
                    dicPAB[locus]=1
                    print locus+" "+str(start_junc)+", "+str(end_junc)+" "+str(read_0)+" "+str(read_3)+" "+str(nb_junc)+" "+str(diclocus[locus][start_junc][end_junc][2])+" "+str(nb_file)
    print "nombre de sn "+str(len(dicsn))
    print "nombre de t "+str(len(dict))
    print "nombre de r "+str(len(dicr))
    print "nombre de new "+str(len(dicnew))
    print "nombre d'autres "+str(len(dicPAB))



