from optparse import OptionParser
from collections import defaultdict

from fastahack import FastaHack

import numpy as np
from sys import stderr
import pandas as pd
import pysam
import time

import cluster
from get_windowed_variance import get_windowed_variance
import genotyper as gt
import IPython
import info_io
import pdb

def get_min_max(cc):
    mn = 9e9
    mx = -1
    for c in cc: 
        s,e = c.get_min_max_start_end()
        mn = min(s,mn)
        mx = max(e,mx)
    
    return mn, mx


if __name__=='__main__':

    opts = OptionParser()

    opts.add_option('','--gglob_dir',dest='gglob_dir')
    opts.add_option('','--dup_tabix',dest='fn_dup_tabix')
    opts.add_option('','--genotype_output',dest='fn_gt_out')
    opts.add_option('','--vcf_output',dest='fn_vcf_out')
    opts.add_option('','--call_output',dest='fn_call_out')
    opts.add_option('','--contig',dest='contig')
    opts.add_option('','--visualizations_dir',dest='out_viz_dir')
    opts.add_option('','--call_table',dest='fn_call_table',default=None)
    opts.add_option('','--total_subsets',dest='total_subsets',type=int,default=1)
    opts.add_option('','--subset',dest='subset',type=int,default=0)
    opts.add_option('','--subset_indivs',dest='subset_indivs', default=None)
    
    opts.add_option('','--genome_fa',dest='fn_fa', default="/net/eichler/vol7/home/psudmant/genomes/fastas/hg19_1kg_phase2_reference/human_g1k_v37.fasta")
    
    opts.add_option('','--do_plot',dest='do_plot',action="store_true",default=False)
    
    opts.add_option('','--simplify_complex_eval',
                    dest='simplify_complex_eval',
                    action="store_true",
                    default=False,
                    help="""complex loci are those where multiple calls overlap
                            The procedure to determine overlappping set of calls and assign
                            each uniquely to a set of individuals is long and works only 
                            really well for high-coverage genomes""")
    
    """
    genotype filtering options
    """
    #opts.add_option('','--singleton_min_sigma',dest='singleton_min_sigma',type=float,default=6.3)
    opts.add_option('','--singleton_min_sigma',dest='singleton_min_sigma',type=float,default=0)
    opts.add_option('','--dup_min_sigma',dest='dup_min_sigma',type=float,default=0) #playing with this didn't work
    opts.add_option('','--filter_min_max_mu_d',dest='min_max_mu_d',type=float,default=0.5)
    opts.add_option('','--filter_max_overlap', dest='max_overlap',type=float,default=0.5)
    opts.add_option('','--max_mu_cp', dest='max_mu_cp',type=float,default=1000)
    opts.add_option('','--filter_X_linked',
                       dest='filter_X_linked',
                       action="store_true",
                       default=False,
                       help="""perform a chi-squared test on calls to see if x-linked and discard if they are""")
    opts.add_option('','--force_output_all', dest='force_output_all', default=False, action="store_true")

    """
    target a distinct ste of regions for genotyping
    """
    opts.add_option('','--target_loci', dest='fn_target_loci',default=None)
    
    (o, args) = opts.parse_args()
    
    subset_indivs = o.subset_indivs
    if subset_indivs!=None:
        subset_indivs = subset_indivs.split(":")
    
    contig = o.contig
    
    target_loci = None
    
    if o.fn_target_loci!=None:
        target_loci = []
        for l in open(o.fn_target_loci):
            c,s,e = l.rstrip().split()
            s,e = int(s), int(e)
            if c == contig: target_loci.append([c,s,e])

    tbx_dups = pysam.Tabixfile(o.fn_dup_tabix)
    callset_clust = cluster.cluster_callsets(o.fn_call_table, contig)
    g = gt.genotyper(contig, gglob_dir=o.gglob_dir, plot_dir=o.out_viz_dir, subset_indivs = subset_indivs, fn_fa=o.fn_fa, dup_tabix = tbx_dups)
    
    fn_sunk_gt_out = o.fn_gt_out.replace(".genotypes",".sunk_genotypes")

    F_gt = open(o.fn_gt_out,'w')
    F_sunk_gt = open(fn_sunk_gt_out,'w')
    F_VCF = open(o.fn_vcf_out,'w')
    F_call = open(o.fn_call_out,'w')
    FINF = open("%s.info"%o.fn_call_out,'w')
    
    info_ob = info_io.info_io(FINF)
    
    g.setup_output(F_gt, F_sunk_gt, F_VCF, F_call, info_ob)
    """
    iterate over lists of overlapping calls
    each element in the list is a recip overlap cluster
    """
    do_plot = o.do_plot
    verbose=False

    k=-1
    
    filt = gt.filter_obj(o.min_max_mu_d, 
                         o.max_overlap,
                         o.max_mu_cp, 
                         o.singleton_min_sigma,
                         o.dup_min_sigma,
                         filter_X_linked = o.filter_X_linked)
     
    for overlapping_call_clusts in callset_clust.get_overlapping_call_clusts(o.total_subsets, o.subset):
        mn, mx = get_min_max(overlapping_call_clusts)
           
        #if contig == "chr2" and not (mx>=75061454 and mn<=75061866): continue
        #if contig == "chr2" and not (mx>=181927132 and mn<=181928581): continue
        #if contig == "chr2" and not (mx>=216225164 and mn<=216226790): continue
        if target_loci:
            for t in target_loci:
                if (mx>=t[1] and mn<=t[2]):
                    break
            else:
                continue

        """
        2 cases
        1. there is one call in the cluster - simply output the call w/ genotypes
        2. the region is more complex
        """
        k+=1

        if k%1==0: 
            print "".join("*" for q in xrange(50))
            print "%d genotypes evaluated..."%k
        
        if len(overlapping_call_clusts) == 1 or o.simplify_complex_eval:
            for c in overlapping_call_clusts:
                start, end = c.get_med_start_end()
                #print start, end
                #if contig == "chr20" and not (start==6648036 and end==6656183): continue
                gt.output(g, contig, start, end, filt, plot=do_plot,v=verbose)  
        else:
            s_e_segs, include_indivs, non_adj = gt.assess_complex_locus(overlapping_call_clusts, g, contig, filt, plot=do_plot)

            
            if len(s_e_segs)<=1 or non_adj:
                for s_e in s_e_segs:
                    s,e = s_e
                    gt.output(g, contig, s, e, filt, plot=do_plot, v=verbose)  
            else:
                for i, s_e in enumerate(s_e_segs):
                    s,e = s_e
                    inc_indivs = include_indivs[i]
                    if o.force_output_all: inc_indivs = None
                    gt.output(g, contig, s, e, filt, include_indivs=inc_indivs, plot=do_plot, v=verbose)  

    F_gt.close()
    F_call.close()
