#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
from tkFileDialog import askopenfilename
from tkFileDialog  import askdirectory
import tkFont
import numpy
import subprocess

class button_browser:
    def __init__(self,master, row, column, name):
        self.file=StringVar()
        self.text=StringVar()
        self.text.set(name)
            #if (name=="Output"):
            #self.text.set("Dossier pour les fichiers de sortie")
        self.label=Label(master, textvariable=self.text).grid(row=row, column=column)
        self.entry=Entry(master, textvariable=self.file, width=30).grid(row=row, column=column+1)
        self.button= Button(master, text="Parcourrir", command = self.browse )
        #x.set(self.browse()))
        self.button.grid(row=row, column=column+3)

    def browse(self):
        Tk().withdraw()
        if (self.text.get()[0]=='D' or self.text.get()=="Output"):
            self.file.set(askdirectory())
        else :
            self.file.set(askopenfilename())

class Window:
    def __init__(self, master):
        y=0
        self.var = IntVar()
        self.multiple = Checkbutton(master, text="Voulez-vous comparer plusieurs fichiers ?",variable=self.var, command = self.test_multiple).grid(row=y, column=1)
        y=y+10
        ### input
        self.input=button_browser(master,y,0,"Fichier d'entree")
        y=y+10
        ### output
        self.output=button_browser(master,y,0,"Prefixe des fichiers de sortie")
        y=y+10
        ### ref
        self.ref=button_browser(master,y,0,"Genome (.fasta)")
        y=y+10
        ### gbk
        self.gbk=button_browser(master,y,0,"Genbank (.gbk)")
        y=y+10
        ### threads
        self.thread_string=StringVar()
        self.thread_string.set("2")
        ref_thread_label=Label(master, text="Nombre de coeurs a utiliser").grid(row=y, column=0)
        ref_thread_entry=Entry(master, textvariable=self.thread_string, width=30).grid(row=y, column=1)
        y=y+10
        self.qbutton= Button(master, text="OK", command=self.Run)
        self.qbutton.grid(row=y, column=3)
    
    #self.entree = Entry(root, textvariable=self.value, width=30).grid(row=20, column=1)
    
    
    def test_multiple(self):
        if (self.var.get()==0):
            self.input.text.set("Fichier d'entree")
            self.input.file.set("")
            self.output.text.set("Prefixe des fichiers de sorties")
            self.output.file.set("")
        else :
            self.input.text.set("Dossier contenant que les fichiers d'entree")
            self.input.file.set("")
            self.output.text.set("Dossier pour les fichiers de sortie")
            self.output.file.set("")


    def Run(self):
        
        if (self.var.get()==0):
            file="run_uniq.sh"
        else:
            file="run_multiple.sh"
        
        fichier=open(file,"w")
        fichier.write("#!/bin/bash\n")
        if (self.var.get()==0):
            fichier.write("INPUT="+self.input.file.get()+"\n")
            fichier.write("OUTPUT="+self.output.file.get()+"\n")
        else:
            fichier.write("INPUTDIR="+self.input.file.get()+"\n")
            fichier.write("OUTPUTDIR="+self.output.file.get()+"\n")
        fichier.write("REF="+self.ref.file.get()+"\n")
        fichier.write("GBK_REF="+self.gbk.file.get()+"\n")
        fichier.write("NB_THREAD="+self.thread_string.get()+"\n")
                          
        fichier.close()
        commande="cat "+file+""
        proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (out, err) = proc.communicate()
        print out
        if (self.var.get()==0):
            commande="cat canvas_run_uniq.sh >> "+file
        else:
            commande="cat canvas_run_multiple.sh >> "+file
        proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (out, err) = proc.communicate()
        commande="bash "+file+ "> test.out"
        proc=subprocess.Popen(commande,shell=True,stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (out, err) = proc.communicate()
        print out
        ###print self.var.get()
        ###print self.input.file.get()
        ###print self.output.file.get()
        ###print self.ref.file.get()
        ###print self.taille_string.get()
        ###print self.gbk.file.get()
        root.quit()


root = Tk()
root.wm_title("Calcul des ARNs circulaires")
window=Window(root)
root.mainloop()


