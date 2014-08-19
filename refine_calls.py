import argparse
import pandas as pd
from sys import stderr
from collections import defaultdict
import math
from scipy import stats
import numpy as np
import cluster

import pysam

from wnd_cp_data import wnd_cp_indiv, dCGH

def output_calls(final_calls, fn):
    """
    output ALL indivdiaul clusters making up calls
    """
    with open(fn, "w") as F:
        #F.write("%s\n"%("\t".join(["indiv_ref", "indiv_test", "chr", "start", "end", "mu",  "p", "window_size"])))
        F.write("%s\n"%("\t".join(["chr", "start", "end", "called_refs, log_likelihood"])))
        for call in final_calls:
            F.write("%s\t%d\t%d\t%s\t%f\n"%(call.contig, call.start, call.end, ",".join(call.ref_indivs),call.total_ll))

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
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--call_table', dest='fn_call_table', required=True, help='Call table file name')
    parser.add_argument('--out_indiv_calls_bed', dest='fn_out_indiv_calls_bed', help='Output filename for individual calls')
    parser.add_argument('--out_clustered_calls_bed', dest='fn_out_clustered_calls_bed')
    parser.add_argument('--out_resolved', dest='fn_out_resolved')
    parser.add_argument('--p_cutoff', type=float, default=0.005, help='P-value call cutoff (Default: %(default)s)')
    parser.add_argument('--min_wnd_call_size', dest='min_wnds', type=int, default=2, help='Min windows included in call (Default: %(default)s)')
    parser.add_argument('--max_callsize', type=int, default=200000, help='Calls must be <%(default)s bp')
    parser.add_argument('--segdups', dest='fn_seg_dups', help='Path to segdups file')
    parser.add_argument('--min_overlapping_calls', type=int, default=3, help='Min number of overlapping dCGH calls (Default: %(default)s)')
    parser.add_argument('--single_window_cutoff', type=float, default=1.0, help='Variance cutoff for single window calls (zero variance messes things up)')
    parser.add_argument('--limit_to_chr', dest='limit_to_chr', default=None, help='List of chrs to include (:-separated, default: %(default)s)')
    parser.add_argument('--indiv_DTS', dest='fn_indiv_DTS', default=None, help='Single sample DTS')
    parser.add_argument('--ref_DTS', dest='fn_ref_DTS', default=None, help='Single reference DTS')
    parser.add_argument('--gglob_dir', default=None, help='Path to gglob directory')
    parser.add_argument('--out_viz_dir', dest='viz_dir', default="./", help='Path to output viz directory (Default: %(default)s)')
    parser.add_argument('--contigs', dest='fn_contigs', default=None, help='Path to tab-delimited table of contigs and contig lengths')
    parser.add_argument('--window_size', type=int, default=None, help='Size of SUNK/wssd sliding windows')
    parser.add_argument('--min_ref_cp_delta', dest='min_d', type=float, default=0, help='Smallest difference in cp between ref and sample to consider (Default: %(default)s)')
    parser.add_argument('--no_P_value_adjust', dest='P_adjust', action='store_false')
    parser.add_argument('--min_mu', default=0.5, type=float, help='Minimum cluster mean (Default: %(default)s)')
    parser.add_argument('--subset_indivs', default=None, help='Colon-separated list of individuals to consider (Default: %(default)s)')

    o = parser.parse_args()    


    subset_indivs = o.subset_indivs
    
    if subset_indivs != None:
        subset_indivs = subset_indivs.split(":")
        subset_indivs = list(set(subset_indivs))

    indiv_DTS = wnd_cp_indiv(o.fn_indiv_DTS, o.fn_contigs, o.window_size) 
    indiv_id = o.fn_indiv_DTS.split("/")[-1].replace("500_bp_","")

    ref_DTSs = {}
    dCGHs = {}
    for fn_ref in o.fn_ref_DTS.split(":"):
        dCGHs[fn_ref.split("/")[-1].replace("500_bp_","")] = dCGH(o.fn_indiv_DTS, 
                                                                  fn_ref,
                                                                  o.fn_contigs,
                                                                  o.window_size)
        
        ref_DTSs[fn_ref.split("/")[-1].replace("500_bp_","")] = wnd_cp_indiv(fn_ref, 
                                                                             o.fn_contigs, 
                                                                             o.window_size) 
        
    call_table = cluster.indiv_callset_table(o.fn_call_table) 
    
    if o.limit_to_chr:
        call_table.filter_by_chr(o.limit_to_chr)
    
    call_table.filter_by_gsize(o.max_callsize)
    call_table.filter(o.p_cutoff, o.min_wnds, o.single_window_cutoff, divide_by_mu=o.P_adjust, min_mu=o.min_mu) 
    call_table.output(o.fn_out_indiv_calls_bed)
    
    tbx_dups = pysam.Tabixfile(o.fn_seg_dups)

    """
    make a call clusterer and use it to get clustered calls
    fn_call_table has all calls made with multiple references
    against the individual
    """
    name = o.fn_call_table.split("/")[-1].split(".")[0]
    call_clusterer = cluster.indiv_cluster_calls(call_table)
    call_clusterer.output_overlap_clusters(o.fn_out_indiv_calls_bed, name)
    final_calls = call_clusterer.resolve_overlapping_clusters(-3, 
                                                              tbx_dups,
                                                              indiv_id,
                                                              indiv_DTS,
                                                              ref_DTSs,
                                                              dCGHs,
                                                              o.gglob_dir,
                                                              o.viz_dir,
                                                              verbose=False, 
                                                              min_overlapping=o.min_overlapping_calls,
                                                              subset_indivs=subset_indivs,
                                                              min_d=o.min_d)

    #output_indiv_clust_elements(final_calls, o.fn_out_resolved)
    output_calls(final_calls, o.fn_out_resolved)
    #output_bed(final_calls, o.fn_out_clustered_calls_bed)
    
