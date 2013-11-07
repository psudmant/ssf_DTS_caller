from optparse import OptionParser
import json
import glob

from sys import stderr
import os
import time
import numpy as np
import pandas as pd

from wnd_cp_data import wnd_cp_indiv
from gglob import gglob

            
if __name__=="__main__":
        
    opts = OptionParser()
    
    opts = OptionParser(usage=usage_txt)
    
    opts.add_option('','--in_DTS_dir',dest='DTS_dir', default=None)
    opts.add_option('','--in_sunk_DTS_dir',dest='sunk_DTS_dir', default=None)

    opts.add_option('','--contigs',dest='fn_contigs', default=None)
    opts.add_option('','--sunk_contigs',dest='fn_sunk_contigs', default=None)
    
    opts.add_option('','--wnd_size',dest='wnd_size', type=int, default=None)
    opts.add_option('','--wnd_slide',dest='wnd_slide', type=int, default=None)
    
    opts.add_option('','--DTS_prefix',dest='DTS_prefix', default="500_bp_")
    
    opts.add_option('','--gglob_dir',dest='gglob_dir')

    (o, args) = opts.parse_args()
    #usage, init, then run
    
    indivs = gglob.get_indivs(DTS_dir, sunk_DTS_dir, DTS_prefix)
    
    contigs = ['chr%d' for d in xrange(1,23)]

    for contig in contigs:   
        
        g = gglob.init_from_DTS(DTS_dir = o.DTS_dir,
                                DTS_prefix = o.DTS_prefix,
                                sunk_DTS_dir = o.sunk_DTS_dir,
                                sunk_DTS_prefix = o.sunk_DTS_prefix,
                                wnd_size = o.wnd_size,
                                indivs = indivs,
                                contig = contig,
                                fn_contigs = o.fn_contigs,
                                fn_sunk_contigs = o.fn_sunk_contigs)
         
