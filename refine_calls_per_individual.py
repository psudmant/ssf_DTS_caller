from optparse import OptionParser
import pandas as pd
from sys import stderr
from collections import defaultdict
import math
from scipy import stats
import numpy as np
import cluster


def output_calls(final_calls, fn):
    
    with open(fn, "w") as F:
        for call in final_calls:
            F.write("%s\n"%call.print_str())

if __name__=="__main__":
    
    opts = OptionParser()
    opts.add_option('', '--call_table', dest='fn_call_table')
    opts.add_option('', '--out_indiv_calls_bed', dest='fn_out_indiv_calls_bed')
    opts.add_option('', '--out_resolved', dest='fn_out_resolved')
    opts.add_option('', '--p_cutoff', dest='p_cutoff', type=float, default=0.001)
    (o, args) = opts .parse_args()
    
    call_table = cluster.callset_table(o.fn_call_table) 
    call_table.filter(o.p_cutoff) 

    """
    make a call clusterer and use it to get clustered calls
    fn_call_table has all calls made with multiple references
    against the individual
    """
    name = o.fn_call_table.split("/")[-1].split(".")[0]
    call_clusterer = cluster.cluster_calls(call_table)
    call_clusterer.output_overlap_clusters(o.fn_out_indiv_calls_bed, name)
    final_calls = call_clusterer.resolve_overlapping_clusters()
    output_calls(final_calls, o.fn_out_resolved)
    

     
