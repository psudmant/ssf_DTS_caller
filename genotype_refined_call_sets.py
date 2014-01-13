from optparse import OptionParser
from collections import defaultdict

import numpy as np
from sys import stderr
import pandas as pd
import pysam
import time

import cluster
from get_windowed_variance import get_windowed_variance
import genotyper as gt




if __name__=='__main__':

    opts = OptionParser()

    opts.add_option('','--gglob_dir',dest='gglob_dir')
    opts.add_option('','--dup_tabix',dest='fn_dup_tabix')
    
    #opts.add_option('','--in_DTS_dir',dest='in_DTS_dir')
    #opts.add_option('','--in_sunk_DTS_dir',dest='in_sunk_DTS_dir')
    #opts.add_option('','--contigs',dest='fn_contigs')
    #opts.add_option('','--sunk_contigs',dest='fn_sunk_contigs')
    #opts.add_option('','--GC_contents_dir',dest='GC_contents_dir')
    opts.add_option('','--genotype_output',dest='fn_gt_out')
    opts.add_option('','--call_output',dest='fn_call_out')
    opts.add_option('','--contig',dest='contig')
    opts.add_option('','--visualizations_dir',dest='out_viz_dir')
    opts.add_option('','--call_table',dest='fn_call_table',default=None)
    opts.add_option('','--total_subsets',dest='total_subsets',type=int,default=1)
    opts.add_option('','--subset',dest='subset',type=int,default=0)
    
    (o, args) = opts.parse_args()
    
    contig = o.contig
    tbx_dups = pysam.Tabixfile(o.fn_dup_tabix)
    callset_clust = cluster.cluster_callsets(o.fn_call_table, contig)
    g = gt.genotyper(contig, gglob_dir=o.gglob_dir, plot_dir=o.out_viz_dir)
    F_gt = open(o.fn_gt_out,'w')
    F_call = open(o.fn_call_out,'w')
    
    """
    iterate over lists of overlapping calls
    each element in the list is a recip overlap cluster
    """
    g.setup_output(F_gt)
    k=-1
    for overlapping_call_clusts in callset_clust.get_overlapping_call_clusts(o.total_subsets, o.subset):
        """
        2 cases
        1. there is one call in the cluster - simply output the call w/ genotypes
        2. the region is more complex
        """
        k+=1

        if k%100==0: print "%d genotypes evaluated..."%k

        if len(overlapping_call_clusts) == 1:
            continue
            c = overlapping_call_clusts[0]
            start, end = c.get_med_start_end()
            gt.output(g, contig, start, end, F_gt, F_call, plot=False)  
        else:
            s_e_segs, include_indivs, non_adj = gt.assess_complex_locus(overlapping_call_clusts, g, contig)
            
            if len(s_e_segs)<=1 or non_adj:
                continue
                for s_e in s_e_segs:
                    s,e = s_e
                    gt.output(g, contig, s, e, F_gt, F_call, plot=False)  
            else:

                for i, s_e in enumerate(s_e_segs):
                    s,e = s_e
                    inc_indivs = include_indivs[i]
                    gt.output(g, contig, s, e, F_gt, F_call, include_indivs=inc_indivs, plot=True)  

    F_gt.close()
    F_call.close()
