from optparse import OptionParser
import json
import glob
import tempfile

from sys import stderr
import sys
import os
import time
import numpy as np
import pandas as pd
import pysam


from wnd_cp_data import wnd_cp_indiv
from gglob import gglob
from wssd_pw_common import *


def fetch_wnds(start, end, wnd_starts, wnd_ends):

    s = np.searchsorted(wnd_starts, start)
    e = np.searchsorted(wnd_ends, end)
    return s,e


class GC_info:

    def __init__(self):
        self.cps = np.array([])
        self.GC = np.array([])

    def add(self, cps, GC):
        self.cps = np.r_[self.cps, cps]
        self.GC = np.r_[self.GC, GC]
    
    def finalize(self, fn_out=None):
        
        if fn_out:
            F = open(fn_out,'w')

        self.cps=np.around(self.cps) 
        d=100.0
        for i in xrange(1,100):
            lower=(i-1)/d
            upper=(i)/d
            
            w = np.where((self.GC>=lower)&(self.GC<upper))
            l = max(w[0].shape[0],1)
            n_cp2 = np.sum(self.cps[w]==2)
            if fn_out:
                print >>F, lower, n_cp2, l, float(n_cp2)/l 
            else:
                print lower, n_cp2, l, float(n_cp2)/l 
    
if __name__=="__main__":
        
    opts = OptionParser()
    
    opts.add_option('','--fn_DTS',dest='fn_DTS', default=None)
    opts.add_option('','--contigs',dest='fn_contigs', default=None)
    opts.add_option('','--wnd_size',dest='wnd_size', type=int, default=None)
    opts.add_option('','--wnd_slide',dest='wnd_slide', type=int, default=None)
    opts.add_option('','--out_dir',dest='out_dir')
    opts.add_option('','--fn_out',dest='fn_out')
    opts.add_option('','--cp2_loci',dest='fn_cp2_loci')
    opts.add_option('','--DTS_prefix',dest='DTS_prefix', default="500_bp_")
    opts.add_option('','--gc_DTS',dest='fn_GC_DTS')

    (o, args) = opts.parse_args()
    #usage, init, then run
    
    GC_DTS = DenseTrackSet(o.fn_contigs,
                           o.fn_GC_DTS,
                           overwrite=False,
                           openMode='r')


    indiv = o.fn_DTS.split("/")[-1].replace("500_bp_","")
    wnd_cp = wnd_cp_indiv(o.fn_DTS, o.fn_contigs, o.wnd_size)
    """
    outstr:
    chr start end indiv 0 0 0 color
    """
    tbx_cp2 = pysam.Tabixfile(o.fn_cp2_loci)
    
    GC_track = GC_info()

    for contig in wnd_cp.contigs:  
        curr_GC_track = GC_info()
        if not contig in tbx_cp2.contigs: continue
        #if not contig in ["chr20"]: continue
        
        print  >>stderr, contig
        
        cps = wnd_cp.get_cps_by_chr(contig)
        wnd_starts, wnd_ends = wnd_cp.get_wnds_by_chr(contig)
        GC =  GC_DTS["GC"][contig][:]
        for tup in tbx_cp2.fetch(contig,parser=pysam.asTuple()):
            cchr, c_s, c_e = tup
            s, e = fetch_wnds(c_s, c_e, wnd_starts, wnd_ends)
            GC_track.add(cps[s:e], GC[s:e])
            curr_GC_track.add(cps[s:e], GC[s:e])
        
        curr_GC_track.finalize()

    GC_track.finalize("%s/%s-%s"%(o.out_dir,o.fn_out,indiv))


