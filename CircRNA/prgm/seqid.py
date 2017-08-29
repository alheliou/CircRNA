#!/usr/bin/env python
# coding=utf8

from Bio import SeqIO
import os
import sys
if __name__ == "__main__":
        rec=SeqIO.read(sys.argv[1],"fasta")
        print rec.id
