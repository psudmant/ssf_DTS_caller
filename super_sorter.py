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
        

class comparator:
    def __init__(self, indivs, group1, group2):
        self.all_indivs = indivs
        
        self.index_by_group = {}
        self.indivs_by_group = {}
        
        
        self.indivs_by_group['group1'] = group1
        self.indivs_by_group['group2'] = group2

        for group, indivs in self.indivs_by_group.iteritems():
            idxs = []
            for indiv in indivs:
                idx = self.all_indivs.index(indiv)
                idxs.append(idx)
            self.index_by_group[group] = np.array(sorted(idxs))
            
            print group, self.index_by_group[group]
    
    def get_coordinates(self, contig, wnd_starts, wnd_ends, poses, mus):

        out_coords = []
        
        if poses.shape[0] == 0: 
            return out_coords
        
        curr_cps = [mus[poses[0]]] 
        out_coords.append([contig, wnd_starts[poses[0]], wnd_ends[poses[0]],np.mean(np.array(curr_cps))])
        for i in xrange(1,poses.shape[0]):
            if poses[i] == poses[i-1]+1: 
                out_coords[-1][2] = wnd_ends[poses[i]]
                curr_cps.append(mus[poses[i]]) 
                out_coords[-1][3] = np.mean(np.array(curr_cps))
            else:
                curr_cps = [mus[poses[i]]] 
                out_coords.append([contig, wnd_starts[poses[i]], wnd_ends[poses[i]], np.mean(np.array(curr_cps))])
        
        return out_coords


    def test_all_dup(self, contig, cps, wnd_starts, wnd_ends, min_cp=3.7):
        """
        test if all indivs in group1 have >cp than all indivs in group2
        """
        mu_all = np.mean(cps[:,:],0)

        poses = np.where(mu_all>min_cp)[0]
        return self.get_coordinates(contig, wnd_starts, wnd_ends, poses, mu_all)
    
    def test_pop_spec_dup(self, contig, cps, wnd_starts, wnd_ends,  group1, group2, min_delta=1):
        """
        test if all indivs in group1 have >cp than all indivs in group2
        """
        idx_1 = self.index_by_group[group1]
        idx_2 = self.index_by_group[group2]

        min_g1 = np.min(cps[idx_1,:],0)
        max_g2 = np.max(cps[idx_2,:],0)
        
        mu_g1 = np.mean(cps[idx_1,:],0)

        delta = min_g1-max_g2
        poses = np.where(delta>min_delta)[0]
        return self.get_coordinates(contig, wnd_starts, wnd_ends, poses, mu_g1)
    
    def test_pop_spec_del(self, contig, cps, wnd_starts, wnd_ends,  group1, group2, min_delta=1):
        """
        test if all indivs in group1 have >cp than all indivs in group2
        """
        idx_1 = self.index_by_group[group1]
        idx_2 = self.index_by_group[group2]

        max_g1 = np.max(cps[idx_1,:],0)
        min_g2 = np.min(cps[idx_2,:],0)
        
        mu_g1 = np.mean(cps[idx_1,:],0)

        delta = min_g2-max_g1
        poses = np.where(delta>min_delta)[0]
        return self.get_coordinates(contig, wnd_starts, wnd_ends, poses, mu_g1)

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

    opts.add_option('','--group1',dest='group1_indivs', default=None)
    opts.add_option('','--group2',dest='group2_indivs', default=None)
    
    opts.add_option('','--outfile',dest='fn_out', default=None)


    opts.add_option('','--test_group1_dup',dest='t_g1_dup', default=False, action='store_true')
    opts.add_option('','--test_group1_del',dest='t_g1_del', default=False, action='store_true')
    opts.add_option('','--test_all_dup',dest='t_all_dup', default=False, action='store_true')
    
    #opts.add_option('','--gglob_dir',dest='gglob_dir')

    (o, args) = opts.parse_args()
    #usage, init, then run
    
    indivs = gglob.get_indivs(o.DTS_dir, o.sunk_DTS_dir, o.DTS_prefix)
    
    indivs1 = o.group1_indivs.split(":")
    indivs2 = o.group2_indivs.split(":")
    if indivs2 == ["-"]: indivs2 = [] 
    comp = comparator(indivs, indivs1, indivs2)
    
    #contigs = ['chr%d'%d for d in xrange(22,23)]
    contigs = ['chr%d'%d for d in xrange(1,23)]
    
    print contigs
    
    #dups is dels too 
    if o.t_g1_dup:
        dups_out = bed_output("%s.dups.wssd.bed"%o.fn_out)
        sunk_dups_out = bed_output("%s.dups.sunk.bed"%o.fn_out)
    
    if o.t_g1_del:
        dels_out = bed_output("%s.dels.wssd.bed"%o.fn_out)
        sunk_dels_out = bed_output("%s.dels.sunk.bed"%o.fn_out)
    
    if o.t_all_dup:
        all_dup_out = bed_output("%s.all_dup.wssd.bed"%o.fn_out)
        all_sunk_dup_out = bed_output("%s.all_dup.sunk.bed"%o.fn_out)
     
    for contig in contigs:   
        #if contig != "chr20": continue
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
        
        if o.t_g1_dup:
            coords = comp.test_pop_spec_dup(contig, g.cp_matrix, g.wnd_starts, g.wnd_ends, "group1", "group2")
            dups_out.add(coords)
            sunk_coords = comp.test_pop_spec_dup(contig, g.sunk_cp_matrix, g.sunk_wnd_starts, g.sunk_wnd_ends, "group1", "group2")
            sunk_dups_out.add(sunk_coords)
        
        if o.t_g1_del:
            coords = comp.test_pop_spec_del(contig, g.cp_matrix, g.wnd_starts, g.wnd_ends, "group1", "group2")
            dels_out.add(coords)
            sunk_coords = comp.test_pop_spec_del(contig, g.sunk_cp_matrix, g.sunk_wnd_starts, g.sunk_wnd_ends, "group1", "group2")
            sunk_dels_out.add(sunk_coords)
        
        if o.t_all_dup:
            coords = comp.test_all_dup(contig, g.cp_matrix, g.wnd_starts, g.wnd_ends)
            all_dup_out.add(coords)
            sunk_coords = comp.test_all_dup(contig, g.sunk_cp_matrix, g.sunk_wnd_starts, g.sunk_wnd_ends)
            all_sunk_dup_out.add(sunk_coords)

    if o.t_g1_dup:
        dups_out.output()
        sunk_dups_out.output()

    if o.t_g1_del:
        dels_out.output()
        sunk_dels_out.output()
    
    if o.t_all_dup:
        all_dup_out.output()
        all_sunk_dup_out.output()

