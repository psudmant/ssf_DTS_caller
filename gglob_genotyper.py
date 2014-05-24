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
        
class bed_output:
    def __init__(self, fn_out):
        self.coords = []
        self.fn_out = fn_out

    def add(self, more_coords):
        self.coords+=more_coords
    
    def output(self):
        with open(self.fn_out,'w') as F:
            for c in self.coords:
                F.write("%s\t%d\t%d\t%f\n"%(c[0],c[1],c[2],c[3]))

def genotype(cps, wnd_starts, wnd_ends, contig, start, end):
    l=cps.shape[1]
    wnd_start = np.searchsorted(wnd_starts, start)   
    wnd_end = min(np.searchsorted(wnd_ends, end)+1,l-1)
    cps = np.mean(cps[:, wnd_start:wnd_end+1],1) 
    return cps

def weighted_genotype(cps, wnd_starts, wnd_ends, contig, start, end):
    """
    Get weighted mean copy numbers for windows that overlap a given region.
    """

    l = cps.shape[1]
    wnd_start = np.searchsorted(wnd_ends[0], start)
    wnd_end = np.searchsorted(wnd_starts[0], end) - 1
    if wnd_start >= l: #Catch if no window overlaps region
        return None
    if wnd_start == wnd_end: #Don't calculate weights if only 1 window overlaps region
        cp_weights = np.array([1])
    else:
    # Calculate % overlap for start and end windows
        start_weight = (max(start, min(end, wnd_ends[0][wnd_start])) \
                      - max(start, wnd_starts[0][wnd_start]))        \
                      / float(wnd_ends[0][wnd_start]-wnd_starts[0][wnd_start])
        end_weight = (min(end, wnd_ends[0][wnd_end]) \
                    - min(end, max(start, wnd_starts[0][wnd_end])))        \
                    / float(wnd_ends[0][wnd_end]-wnd_starts[0][wnd_end])
        cp_weights = np.array([start_weight] + [1 for x in xrange(wnd_start + 1, wnd_end)] +
                              [end_weight])
    copies = np.average(cps.ix[:, wnd_start:wnd_end], axis = 1, weights = cp_weights)
    return copies

if __name__=="__main__":
         
    opts = OptionParser()
    
    opts.add_option('','--in_DTS_dir',dest='DTS_dir', default=None)
    opts.add_option('','--in_sunk_DTS_dir',dest='sunk_DTS_dir', default=None)

    opts.add_option('','--contigs',dest='fn_contigs', default=None)
    opts.add_option('','--sunk_contigs',dest='fn_sunk_contigs', default=None)
    
    opts.add_option('','--wnd_size',dest='wnd_size', type=int, default=None)
    opts.add_option('','--wnd_slide',dest='wnd_slide', type=int, default=None)
    
    opts.add_option('','--DTS_prefix',dest='DTS_prefix', default="500_bp_")
    opts.add_option('','--sunk_DTS_prefix',dest='sunk_DTS_prefix', default="500_bp_")
    
    #opts.add_option('','--inf',dest='fn_inf', default=None)
    #opts.add_option('','--outdir',dest='outdir', default=None)

    opts.add_option('','--outfile',dest='fn_out', default=None)
    opts.add_option('','--input_loci',dest='fn_input_loci', default=None)


    #opts.add_option('','--gglob_dir',dest='gglob_dir')

    (o, args) = opts.parse_args()
    #usage, init, then run
    
    indivs = gglob.get_indivs(o.DTS_dir, o.sunk_DTS_dir, o.DTS_prefix)
    
    t_loci = pd.read_csv(o.fn_input_loci, header=None, names=["chr", "start", "end"], delimiter="\t")
    t_loci = t_loci.sort(["chr", "start", "end"])
    
    
    F_cps = open("%s.genotypes"%o.fn_out,'w')
    F_sunk_cps = open("%s.sunk_genotypes"%o.fn_out,'w')
    
    F_cps.write("%s\n"%("\t".join(["locus"] + indivs)))
    F_sunk_cps.write("%s\n"%("\t".join(["locus"] + indivs)))
    
    for contig, t_by_chr in t_loci.groupby('chr', sort=False): 
        #if contig !="chr20": continue
        print >>stderr, contig
        
        g = gglob.init_from_DTS(DTS_dir = o.DTS_dir,
                                DTS_prefix = o.DTS_prefix,
                                sunk_DTS_dir = o.sunk_DTS_dir,
                                sunk_DTS_prefix = o.sunk_DTS_prefix,
                                wnd_size = o.wnd_size,
                                wnd_slide = o.wnd_slide,
                                indivs = indivs,
                                contig = contig,
                                fn_contigs = o.fn_contigs,
                                fn_sunk_contigs = o.fn_sunk_contigs)
        
        cps_by_locus = {}
        sunk_cps_by_locus = {}
        loci = [] 
        for idx, row in t_by_chr.iterrows(): 
            print row['chr'], row['start'], row['end'] 
            cps = genotype(g.cp_matrix, g.wnd_starts, g.wnd_ends, row['chr'], row['start'], row['end']) 
            sunk_cps = genotype(g.sunk_cp_matrix, g.sunk_wnd_starts, g.sunk_wnd_ends, row['chr'], row['start'], row['end']) 
            key = "%s:%d-%s"%(row['chr'], row['start'], row['end']) 
            cps_by_locus[key] = cps
            sunk_cps_by_locus[key] = sunk_cps
            loci.append(key)
        
        for locus in loci:
           cps = cps_by_locus[locus]
           F_cps.write("%s\t%s\n"%(locus,"\t".join(["%f"%cp for cp in cps])))
           sunk_cps = sunk_cps_by_locus[locus]
           F_sunk_cps.write("%s\t%s\n"%(locus,"\t".join(["%f"%cp for cp in sunk_cps])))

