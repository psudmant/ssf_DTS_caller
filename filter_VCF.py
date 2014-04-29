from optparse import OptionParser
import pandas as pd
import pysam
from sys import stderr
from sys import stdin

from csv import DictReader
import gzip
import subprocess


import cStringIO
io_method = cStringIO.StringIO


if __name__=="__main__":

    opts = OptionParser()
    opts.add_option('','--DEL_lscore', dest='del_lscore', type = float)
    opts.add_option('','--DUP_lscore', dest='dup_lscore', type = float)
    opts.add_option('','--CNV_lscore', dest='CNV_lscore', type = float)
    o, args = opts.parse_args()
    
    type_to_cutoff = {"DEL":o.del_lscore,
                     "DUP":o.dup_lscore,
                     "CNV":o.CNV_lscore}
    for l in stdin:
        if l[0:1] == "#":
            print l.rstrip()
        else:
            t = l.split("\t",8)

            loc = t[2] 
            inf =  t[7]
            inf_d = {i.split('=')[0]:i.split("=")[1]  for i in inf.split(";") }
            if float(inf_d["LPROBS"]) >= type_to_cutoff[inf_d['SVTYPE']]:
                print l.rstrip()

