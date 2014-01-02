from optparse import OptionParser
from collections import defaultdict

import numpy as np
from sys import stderr
import pandas as pd


from get_windowed_variance import get_windowed_variance
from genotyper import genotyper
from genotyper import genotype_table 



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
    opts.add_option('','--call_table',dest='fn_call_table',default=None)
    
    (o, args) = opts.parse_args()
        
    tbx_dups = pysam.Tabixfile(o.fn_dups)
    
    #load up the calls
    #sort and 




    for contig in tbx_gt_loci.contigs:
        print >>stderr, "working on contig:%s"%contig
        if contig!="chr10": continue
         
        #g = genotyper(o.in_DTS_dir, o.in_sunk_DTS_dir, o.fn_contigs, o.fn_sunk_contigs, o.wnd_size, contig, o.out_viz_dir, limit_to_n=1000)
        #g = genotyper(o.in_DTS_dir, o.in_sunk_DTS_dir, o.fn_contigs, o.fn_sunk_contigs, o.wnd_size, contig, o.out_viz_dir, limit_to_n=20)
        g = genotyper(contig, gglob_dir=o.gglob_dir, plot_dir=o.out_viz_dir)
        FOUT = open("%s/%s.txt"%(o.out_viz_dir,contig),'w')
                      
        #g = genotyper(o.in_DTS_dir, None, o.fn_contigs, None, o.wnd_size, contig, o.out_viz_dir, limit_to_n=400)
        g_table = genotype_table(g)

        for locus in tbx_gt_loci.fetch(contig,parser=pysam.asTuple()):
            chr, start, end, p_val = locus
            start, end = int(start), int(end)
            #if start != 56445975: continue
            print chr, start, end 
            X = g.get_gt_matrix(chr, start, end)
            Xs = g.get_sunk_gt_matrix(chr, start, end)
            FOUT.write("%s %d %d "%(chr,start,end)) 
            gmmX, labelsX, Z, vbnds_X  = g.GMM_genotype(X,FOUT)
            FOUT.write(" ")
            gmmXs, labelsXs, Zs, vbnds_Xs  = g.GMM_genotype(Xs,FOUT)
            FOUT.write("\n ")
            g.plot(X, Xs, gmmX, gmmXs, vbnds_X, vbnds_Xs, Z, Zs, chr, start, end)
            #g.plot(X, Xs, gmmX, gmmXs, vbnds_X, vbnds_Xs, chr, start, end)
             
            #labels, n_clusts, s_score  = g.ms_genotype(X)
            #labels, n_clusts, s_score  = g.AP_genotype(X)
            #labels, n_clusts, s_score  = g.DBS_genotype(X)
            #if n_clusts == -1 or n_clusts == 1: continue
            #g.plot(X, labels, chr, start, end)
            
            #g_table.add(chr, start, end, X, labels, s_score)
        g_table.output("%s/%s.vcf"%(o.outdir,contig))     

    exit(1)
    #cp_data = dCGH(o.fn_in_DTS,o.fn_ref_DTS,o.fn_contigs,wnd)    
    
    null_dist = null_distribution(tbx_gc)

    segment_callset = callset()
    
    caller_by_chr = {}
    
    for chr in chrs:
        print >>stderr,"%s..."%chr
        magnitude_vect = cp_data.get_cps_by_chr(chr)
        starts_vect,ends_vect = cp_data.get_wnds_by_chr(chr) 

        #plot_GC(chr,tbx_gc,magnitude_vect,starts_vect,ends_vect)
        
        gapped_wnds = cp_data.get_overlapping_wnds(chr,tbx_gaps) 
        segdup_wnds = cp_data.get_overlapping_wnds(chr,tbx_dups) 
        null_dist.add(magnitude_vect,[gapped_wnds,segdup_wnds])
        
        print>>stderr, "cutoff_scale:%f"%( cutoff_scale )
        caller = ssf_caller(chr,magnitude_vect,starts_vect,ends_vect,cutoff_scale,use_means=True,max_merge=max_merge,scale_width=1)
        caller_by_chr[chr] = caller
        
        #plotter = line_plot(chr,caller,o.fn_gene_tabix,o.fn_dup_tabix,o.plot_lims,o.out_viz_dir) 
        #segment_callset_chr[chr] = 
        segment_callset += caller.get_callset([tbx_gaps]) 
    
    segment_callset.get_p_values(null_dist)
    #segment_callset.get_t_stats(null_dist)
    #segment_callset.get_bh_corrected_significant(o.fdr)
    
    if o.do_plot: 
        for chr, caller in caller_by_chr.iteritems():
            plotter = plot(chr,caller,o.fn_gene_tabix,o.fn_dup_tabix,o.plot_lims,o.out_viz_dir,segment_callset) 
            #plotter.plot_all(chunk_len=1000000,bp_start=150000000,add_heatmap=False)
            plotter.plot_all(chunk_len=1000000,bp_start=0,add_heatmap=False)

    fn_out = "%s/%s.segs.bed"%(o.out_call_dir,o.out_prefix)
    segment_callset.output(fn_out,t_stats=False)

    #null_dist.explore_mu_based_p_values()
    #for chr in chrs:
    #    for seg in segments_by_chr[chr]:
    #        start,end,val,wnd_start,wnd_end = seg
    #        wnd_width=wnd_end-wnd_start+1
    #        p_val = nd.get_p_value(wnd_width,val)
    #        print >>OUTF, "%s\t%d\t%d\t%f\t%f"%(chr,seg[0],seg[1],seg[2],p_val)



