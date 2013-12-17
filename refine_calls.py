from optparse import OptionParser
import pandas as pd
from sys import stderr
from collections import defaultdict
import math
from scipy import stats
import numpy as np
import cluster

import pysam


def output_calls(final_calls, fn):
    """
    output ALL indivdiaul clusters making up calls
    """
    with open(fn, "w") as F:
        F.write("%s\n"%("\t".join(["indiv_ref", "indiv_test", "chr", "start", "end", "mu",  "p", "window_size"])))
        for call in final_calls:
            F.write("%s"%call.print_str())

def output_bed(final_calls, fn):
    indiv  = fn.split("/")[-1].split(".")[0]
    with open(fn, "w") as F:
        F.write("""track name="%s_merged" description="%s_merged" visibility=2 itemRgb="On"\n"""%(indiv, indiv))
        for call in final_calls:
            F.write("%s"%call.bed_str())

def output_indiv_clust_elements(final_calls, fn):
    """
    output ALL indivdiaul calls
    inside clusters to be combined at a later
    date
    """
    with open(fn, "w") as F:
        F.write("%s\n"%("\t".join(["indiv_ref", "indiv_test", "chr", "start", "end", "mu",  "p", "window_size"])))
        for call in final_calls:
            F.write("%s"%call.get_indiv_calls_str())

if __name__=="__main__":
    
    opts = OptionParser()
    opts.add_option('', '--call_table', dest='fn_call_table')
    opts.add_option('', '--out_indiv_calls_bed', dest='fn_out_indiv_calls_bed')
    opts.add_option('', '--out_clustered_calls_bed', dest='fn_out_clustered_calls_bed')
    opts.add_option('', '--out_resolved', dest='fn_out_resolved')
    opts.add_option('', '--p_cutoff', dest='p_cutoff', type=float, default=0.005)
    opts.add_option('', '--min_wnd_call_size', dest='min_wnds', type=int, default=2)
    opts.add_option('', '--max_callsize', dest='max_callsize', type=int, default=200000)
    opts.add_option('', '--segdups', dest='fn_seg_dups')
    opts.add_option('', '--min_overlapping_calls', dest='min_overlapping_calls', type=int, default=1)
    opts.add_option('', '--single_window_cutoff', dest='single_window_cutoff', type=float, default=1.0)
    """
        min_overlapping_calls in the minimum # of calls

    """
    (o, args) = opts .parse_args()
    
    call_table = cluster.callset_table(o.fn_call_table) 
    call_table.filter_by_gsize(o.max_callsize)
    call_table.filter(o.p_cutoff, o.min_wnds,single_window_cutoff, divide_by_mu=True) 
    call_table.output(o.fn_out_indiv_calls_bed)
    
    tbx_dups = pysam.Tabixfile(o.fn_seg_dups)

    """
    make a call clusterer and use it to get clustered calls
    fn_call_table has all calls made with multiple references
    against the individual
    """
    name = o.fn_call_table.split("/")[-1].split(".")[0]
    call_clusterer = cluster.cluster_calls(call_table)
    #call_clusterer.output_overlap_clusters(o.fn_out_indiv_calls_bed, name)
    final_calls = call_clusterer.resolve_overlapping_clusters(-6, tbx_dups, verbose=False, min_overlapping=o.min_overlapping_calls)
    #output_indiv_clust_elements(final_calls, o.fn_out_resolved)
    output_calls(final_calls, o.fn_out_resolved)
    #output_bed(final_calls, o.fn_out_clustered_calls_bed)
    

     
