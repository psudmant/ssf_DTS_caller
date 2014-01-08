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
    
    """
        HERE, these shoudl be recip overlap too
        ALSO, you may need to refine them?
    """
    for c in callset_clust.get_overlapping_call_clusts(o.total_subsets, o.subset):
        gt.get_best_gt(c, contig, g)
