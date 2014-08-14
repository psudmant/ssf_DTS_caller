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

    opts.add_option('','--fn_gts',dest='fn_gts')
    opts.add_option('','--fn_out',dest='fn_out')
    opts.add_option('','--exclude_indivs',dest='exclude_indivs', default=None)
    
    (o, args) = opts.parse_args()
    
    if o.exclude_indivs is not None:
        exclude_indivs = o.exclude_indivs.split(":")
    else:
        exclude_indivs = []
    
    FIN = open(o.fn_gts)
    FOUT = open(o.fn_out,'w')
    

    h = FIN.readline()
    h_indivs = h.rstrip().split("\t")[3:]
    
    new_inds = []
    new_inds_idxs = []

    for i, ind in enumerate(h_indivs):
        if not ind in exclude_indivs:
            new_inds.append(ind)
            new_inds_idxs.append(i)
    
    indivs = new_inds
    idxs = np.array(new_inds_idxs)
    FOUT.write("contig\tstart\tend\t{inds}\n".format(inds="\t".join(indivs)))
    

    af_counts = defaultdict(int)
    af_gt50k = defaultdict(int)

    for l in FIN:
        gts = np.array([int(g) for g in l.rstrip().split()[3:]])
        c,start,end = l.split("\t")[0:3]
        gts = gts[idxs]
        s = 0
        
        if np.sum(gts==2) == gts.shape[0]:continue

        
        for i in xrange(3):
            s+=np.sum(gts==i)
        
        if s==gts.shape[0] and np.sum(gts==gts[0])!=gts.shape[0]:
            FOUT.write("{chr}\t{s}\t{e}\t{gts}\n".format(chr=c,
                                                         s=start,
                                                         e=end,
                                                         gts="\t".join(["%d"%g for g in gts])))
            ac = (2*np.sum(gts==0)+np.sum(gts==1))
            af = float(ac)/(2*gts.shape[0])
            af = af>0.5 and 1-af or af
            af_counts[af]+=1
            if int(end)-int(start)>50000:
                af_gt50k[af]+=1
            
    FOUT.close()
    FIN.close()
    

    FOUT = open("%s.afs"%o.fn_out,'w')
    print >> FOUT, "allele frequency\tcount\tsize_range"
    for af, count in af_counts.iteritems():
        print >> FOUT, "%f\t%d\tall"%(af,count)
        print >> FOUT, "%f\t%d\tgt_50kb"%(af,af_gt50k[af])

    FOUT.close()
