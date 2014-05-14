import glob
from sys import stderr
from sys import stdout
import numpy as np

import time
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mCols
import matplotlib.cm as cm
import matplotlib.mlab as mlab

##local repo
from sklearn import cluster 
from sklearn import metrics
from sklearn import mixture

from fastahack import FastaHack

from sets import Set

import math
import random
from scipy.stats import norm
from scipy.stats import fisher_exact

import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

from wnd_cp_data import wnd_cp_indiv
from gglob import gglob

import cluster as m_cluster
from sets import Set
from scipy.stats.mstats import mode

import IPython
import pdb


class call:

    def __init__(self, chr, start, end, X, labels, s):
        self.chr = chr
        self.start = start
        self.end = end
        self.X = X
        self.labels = labels
        self.s = s

class genotype_table:

    def __init__(self, gtyper):
        self.gtyper = gtyper
        self.calls = []
    
    def add(self, chr, start, end, X, labels, s):
        self.calls.append(call(chr, start, end, X, labels, s))
    
    def output(self, fn):
        
        indiv_str = "\t".join(self.gtyper.indivs)
        header="chr\tstart\tend\ts\t%s"%indiv_str
        with open(fn, 'w') as F:
            for call in self.calls:
                F.write("%s\t%d\t%d\t%f\t%s\n"%(call.chr,
                                call.start,
                                call.end,
                                call.s,
                                "\t".join([str(l) for l in call.labels])))




def get_correlation_matrix(starts_ends, g, contig, outdir=None, do_plot=False, append=""):
    
    n_indivs = len(g.indivs)
    l = len(starts_ends)-1
    mus = np.zeros((n_indivs,l))

    starts_ends = sorted(np.unique(starts_ends))

    n_indivs = len(g.indivs)
    l = len(starts_ends)-1
    mus = np.zeros((l,n_indivs))
    
    min_s, max_e = starts_ends[0], starts_ends[-1]
    for i in xrange(l):
        s, e = starts_ends[i], starts_ends[i+1] 
        X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        mus[i,:] = np.mean(X,1)   
        print contig, s, e
    
    c =  np.corrcoef(mus)
    
    #zero out diagnals
    off_diag = []
    for j in xrange(c.shape[0]-1):
        off_diag.append(c[j,j+1])
    
    if do_plot:
        print "PLAAAAAAATING!"
        plt.gcf().set_size_inches(14,6)
        fig, axes = plt.subplots(2)
        p = axes[0].pcolor(c)
        fig.colorbar(p, cax=axes[1])
        plt.savefig("%s/%s_cor_%d_%d_bps%s.png"%(outdir, contig, min_s, max_e, append))
        """
        #PLOT INDIVIDUAL CORRELATIONS
        plt.gcf().clear()
        
        d = int(np.floor(np.sqrt(l)))+1
        plt.gcf().set_size_inches(18,10)
        
        fig, axes = plt.subplots(d,d)
        
        font = {'family' : 'normal', 'weight': 'normal', 'size': 4} 
        matplotlib.rc('font', **font)
        
        for i in xrange(d): 
            for j in xrange(d): 
                mus[0]
                #p = (j*d)+i
                p = (i*d)+j
                if p<l-1:
                    axes[i,j].plot(mus[p], mus[p+1],'b.', ms=.1)
                    min_x, max_y = np.amin(mus[p]), (np.amin(mus[p+1])+np.amax(mus[p+1]))/2.0
                    axes[i,j].text(min_x,max_y,"%.2f"%off_diag[p], fontsize=4)
                    axes[i,j].text(min_x,max_y-.5,"%d"%starts_ends[p+1], fontsize=4)
                    axes[i,j].xaxis.set_tick_params(width=1)
                    axes[i,j].yaxis.set_tick_params(width=1)
        plt.savefig("%s/%s_i_cors_%d_%d_bps.pdf"%(outdir, contig, min_s, max_e))
        plt.gcf().clear()
        """

    return c, off_diag, mus

def get_correlated_segments(all_starts_ends, g, contig, f_r_cutoff, outdir, do_plot=False):
    """
    take a set of coordinates representing starts and ends
    over a locus and cluster them into contiguous chunks that
    likely represent the same underlying call as a result of their
    high correlation
    """
    
    all_starts_ends = sorted(np.unique(all_starts_ends))
    
    c, off_diag, mus = get_correlation_matrix(all_starts_ends, g, contig, outdir=outdir, do_plot=do_plot)
    
    #merge based on a sliding scale based on copy 
    r_cutoff = f_r_cutoff(np.mean(mus)) 
    #print all_starts_ends 
    original_c = c
    prev_gts = None
    count = 0
    while(np.amax(off_diag)>r_cutoff):
        count+=1
        to_pop = []
        for i in xrange(len(all_starts_ends)-2):
            if np.absolute(off_diag[i])>=r_cutoff:
                to_pop.append(i+1)
        new_positions = [] 
        for i, v in enumerate(all_starts_ends):
            if not i in to_pop:
                new_positions.append(v)
       
        all_starts_ends = np.unique(new_positions)
        #print off_diag
        #print all_starts_ends 
        if len(all_starts_ends) == 2: break
        c, off_diag, mus = get_correlation_matrix(all_starts_ends, g, contig, outdir=outdir, do_plot=do_plot, append = "_merge_%d"%count)

    s_e_tups = []
    for i in xrange(len(all_starts_ends)-1):
        s_e_tups.append([all_starts_ends[i], all_starts_ends[i+1]])
    
    pdb.set_trace()
    return s_e_tups, original_c, mus

def cluster_segs(segs, max_frac_uniq=0.2):
    clust = m_cluster.call_cluster()
    for seg in segs:
        s,e = seg
        clust.add({'start':s, 'end':e})
        
    groups = clust.cluster_by_recip_overlap(0.85, max_frac_uniq=max_frac_uniq)
    
    merged_segs = [] 
    for g in groups:
        sub_segs = [[c['start'], c['end']] for c in g.calls]
        merged_segs.append(sub_segs) 
    
    return merged_segs
  

def get_uniq_chunks(segs, indivs_by_cnv_segs, indivs):
    
    s_indivs = Set(indivs)
    segs = sorted(segs, key=lambda x:(x[0],x[1]))
    mn,mx = segs[0][0],segs[0][1]
    uniq_chunks = [[]] 
    indivs_in_chunk = [[]]
    
    curr_indivs = [] 
    for seg in segs:
        if seg[0]<=mx and seg[1]>=mn:
            uniq_chunks[-1].append(seg)
            indivs_in_chunk[-1] += indivs_by_cnv_segs[tuple(seg)]
            mn,mx = min(mn,seg[0]), max(seg[1],mx)
        else:
            indivs_in_chunk.append([])
            uniq_chunks.append([])
            uniq_chunks[-1].append(seg)
            indivs_in_chunk[-1] += indivs_by_cnv_segs[tuple(seg)]
            mn,mx = seg
    
    n_indivs_in_chunk = []
    for inds in indivs_in_chunk:
        #take indivs AND indivs in the chunk
        n_indivs_in_chunk.append(list(Set(inds).intersection(s_indivs)))

    return uniq_chunks, n_indivs_in_chunk

def get_segs_from_clustered_indivs(indivs_by_grp, segs_by_grp, indivs_by_cnv_segs, g):
    
    indivs_by_segs = {}

    for grp, indivs in indivs_by_grp.iteritems():
        segs = segs_by_grp[grp]
        #print len(indivs), segs
        uniq_chunks,indivs_in_chunks = get_uniq_chunks(segs, indivs_by_cnv_segs,indivs)
        for i, u_chunk_ses in enumerate(uniq_chunks):
            #print '\t', u_chunk_ses, indivs
            med_s = int(np.median(np.array([se[0] for se in u_chunk_ses])))
            med_e = int(np.median(np.array([se[1] for se in u_chunk_ses])))
            if len(indivs_in_chunks[i]) == len(indivs):
                indivs_by_segs[tuple([med_s,med_e])] = indivs
            #print "\t", med_s, med_e, len(indivs_in_chunks[i])
         
    return indivs_by_segs

def cluster_indivs(cnv_segs_by_indiv, g, cutoff=0.85, max_frac_uniq=0.1):
    
    l = len(g.indivs)
    mat = np.zeros([l,l])

    for i, indiv_1 in enumerate(g.indivs):
        for j, indiv_2 in enumerate(g.indivs):
            segs1 = cnv_segs_by_indiv[indiv_1] 
            segs2 = cnv_segs_by_indiv[indiv_2]
             
            if len(segs1) == 0 or len(segs2) ==0:
                mat[i,j] = 1
            elif i==j:
                mat[i,j] = 0 
            else:
                mat[i,j] = 1-m_cluster.seg_sets_intersect(cnv_segs_by_indiv[indiv_1], cnv_segs_by_indiv[indiv_2], max_frac_uniq)
    
    grps = m_cluster.linkage_cluster(mat, cutoff)
    
    indivs_by_clust = {}
    segs_by_clust = {}
    
    for grp in np.unique(grps):
        inds = []
        segs = []
        for j in np.where(grps==grp)[0]:
            if len(cnv_segs_by_indiv[g.indivs[j]])>0:
                inds.append(g.indivs[j])
                segs+=cnv_segs_by_indiv[g.indivs[j]]
        if len(inds)>0: 
            indivs_by_clust[grp] = inds
            segs_by_clust[grp] = segs

    return indivs_by_clust, segs_by_clust 

def clust_seg_groups_by_gt(clust, indivs, g, indivs_by_cnv_segs, contig):
    
    all_gts = []
    for seg in clust:
        s,e = seg
        X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        gX = g.GMM_genotype(X)
        #g_by_indiv = gX.get_gts_by_indiv()
        g_by_indiv, gts_to_label, label_to_gt = gX.get_gts_by_indiv()

        gts = [g_by_indiv[indiv] for indiv in indivs]
        all_gts.append(np.array(gts))
    
    gt_clusts = [{"gts":all_gts[0], "segs":[clust[0]]}]
    
    for i in xrange(1,len(clust)):
        for gt_clust in gt_clusts:
            if np.all(gt_clust["gts"]==all_gts[i]):
                gt_clust["segs"].append(clust[i])
                break
        else:          
            gt_clusts.append({"gts":all_gts[i], "segs":[clust[i]]})
    
    final_seg_to_indiv = {}
    for gt_clust in gt_clusts: 
        max_indivs = -1
        best_seg = []
        all_indivs = []
        for seg in gt_clust['segs']:
            if len(indivs_by_cnv_segs[tuple(seg)])>max_indivs:
                max_indivs=len(indivs_by_cnv_segs[tuple(seg)])
                best_seg = tuple(seg)
                all_indivs+=indivs_by_cnv_segs[tuple(seg)]
            
        final_seg_to_indiv[best_seg] = list(Set(all_indivs)) 
        
    return final_seg_to_indiv

def cluster_overlapping_idGTs(indivs_by_cnv_segs, g, contig, max_uniq_thresh):

    segs = []
    for seg, indivs in indivs_by_cnv_segs.iteritems():
        segs.append(seg)
    
    clusts = cluster_segs(segs, max_uniq_thresh)
    
    new_inds_by_seg = {}

    for clust in clusts:
        if len(clust)==1:
            tup = tuple(clust[0])
            new_inds_by_seg[tup] = indivs_by_cnv_segs[tup]
        
        inds_by_seg = clust_seg_groups_by_gt(clust, g.indivs, g, indivs_by_cnv_segs, contig)
        
        for seg, inds in inds_by_seg.iteritems():
            new_inds_by_seg[tuple(seg)] = inds
            
    return new_inds_by_seg
                
     
def assess_complex_locus(overlapping_call_clusts, g, contig, filt, r_cutoff = lambda x: (0.9/x)+0.5, plot=False):
    """
    First chop up into ALL constituate parts
    """
    
    all_starts_ends = []
    min_s, max_e = 9e9, -1
    for clust in overlapping_call_clusts: 
        s,e = clust.get_med_start_end()
        min_s = min(min_s, s)
        max_e = max(max_e, e)
        all_starts_ends.append(s)
        all_starts_ends.append(e)
    
    """
    #merge correllated calls
    #commented for now...
    """

    s_e_segs, c, mus = get_correlated_segments(all_starts_ends, g, contig, r_cutoff, "./plotting/test", do_plot=plot)
   
    mu_cp = np.mean(mus) 
    
    print "%d segs merged to %d correlated segs"%(len(all_starts_ends), len(s_e_segs))
    #instead of correlation cleaning, use below 4
    #all_starts_ends = sorted(np.unique(all_starts_ends))
    #s_e_segs = []
    #for i in xrange(len(all_starts_ends)-1):
    #    s_e_segs.append([all_starts_ends[i], all_starts_ends[i+1]])

    t = time.time()
    """
    now remove the non-varying chunks
    """
    print "getting var chunks..."
    #THIS IS EATING UP TIME
    CNV_segs, CNV_gXs = filter_invariant_segs(s_e_segs, g, contig) 
    print "got var chunks in %fs"%(time.time()-t)
    
    if len(CNV_segs) <= 1 or non_adjacent(CNV_segs):
        indivs_to_assess = [None for i in s_e_segs]
        exclude_loci = [None for i in s_e_segs]
        return s_e_segs, indivs_to_assess, True
    
    else:

        cnv_segs_by_indiv = {}
        for i, indiv in enumerate(g.indivs):
            indiv_cnv_segs = []
            for i, gX in enumerate(CNV_gXs):
                seg = list(CNV_segs[i])
                if gX.is_var(indiv, g, force_not_mode=True):
                    if len(indiv_cnv_segs)>0 and seg[0] == indiv_cnv_segs[-1][1]: 
                        indiv_cnv_segs[-1][1] = seg[1]
                    else:
                        indiv_cnv_segs.append(seg)
            cnv_segs_by_indiv[indiv] = indiv_cnv_segs
        
        indivs_by_cnv_segs = {}
        for indiv, cnv_segs in  cnv_segs_by_indiv.iteritems():
            for s_e in cnv_segs:
                s_e_tup = tuple(s_e)
                if not s_e_tup in indivs_by_cnv_segs:
                    indivs_by_cnv_segs[s_e_tup] = []
                indivs_by_cnv_segs[s_e_tup].append(indiv)
        
        """
        if high copy, then force more clustering

        max_uniq_thresh: this is the cluster threshold for the clustering of individual segments together

        max_frac_uniq: when clustering individuals by their similar sets of overlapping calls, this is 
        the fraction of below which if two individuals are different we just group them together
        this has to happen for high copy regions because they are so nasty
        """
        
        max_uniq_thresh=0.3
        max_frac_uniq=.1

        if mu_cp>10:
            max_frac_uniq=.5
            max_uniq_thresh=0.7
        elif mu_cp>20:
            max_frac_uniq=.9
            max_uniq_thresh=0.85

        indivs_by_grp, segs_by_grp = cluster_indivs(cnv_segs_by_indiv,g, 0.85, max_frac_uniq=max_frac_uniq) 
        indivs_by_cnv_segs = get_segs_from_clustered_indivs(indivs_by_grp, segs_by_grp, indivs_by_cnv_segs, g)
        
        #finally cluster these

        if len(indivs_by_cnv_segs.keys())>1:
            indivs_by_cnv_segs = cluster_overlapping_idGTs(indivs_by_cnv_segs, g, contig, max_uniq_thresh)
        
        if plot: 
            m_cluster.cluster_callsets.plot(overlapping_call_clusts, "./plotting/test", g, indivs_by_cnv_segs, [], CNV_segs, cnv_segs_by_indiv)
        
        #pdb.set_trace()
        #m_cluster.cluster_callsets.plot(overlapping_call_clusts, "./plotting/test/", g, c, s_e_segs, CNV_segs, cnv_segs_by_indiv) 
        s_e_segs = []
        indivs_to_assess = []
        
        for cnv_seg, inds_called in indivs_by_cnv_segs.iteritems():
            assess_indivs = Set(g.indivs)   
            s1, e1 = cnv_seg
            for cnv_seg2, indivs2 in indivs_by_cnv_segs.iteritems():
                s2, e2 = cnv_seg2
                if cnv_seg != cnv_seg2 and overlap(s1, e1, s2, e2):
                    assess_indivs = assess_indivs - Set(indivs2)
            
            s_e_segs.append(cnv_seg) 
            indivs_to_assess.append(list(assess_indivs))
        
        return s_e_segs, indivs_to_assess, False
   
def filter_invariant_segs(s_e_segs, g, contig):

    """
       remove crappy segs
       filter the segs to be only regions that are CNV
       AND
       merge adjacent calls with equal genotypes (technically the above should handle
       this, but, occationally it doesn't)
    """

    CNV_segs = []
    CNV_gXs = []
    prev_labels = None
    for s_e in s_e_segs:
        s, e = s_e
        X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        gX = g.GMM_genotype(X)
        #print s, e, gX.n_clusts
        if gX.n_clusts > 1:
            if np.all(prev_labels == gX.labels):
                CNV_segs[-1][1] = e 
            else:
                CNV_segs.append([s,e])
                CNV_gXs.append(gX)
          
        prev_labels = gX.labels
    
    return CNV_segs, CNV_gXs        

    
def overlap(s1, e1, s2, e2):
    if (s1<=e2 and s2<=e1):
        return True
    return False

def non_adjacent(CNV_segs):
    
    for i in xrange(len(CNV_segs)-1):
        if CNV_segs[i][1] == CNV_segs[i+1][0]:
            return False

    return True

def test_correlation(overlapping_call_clusts, g, contig):
    
    ##ONLY LOOK AT INDIVS w/ calls???
    all_starts_ends = []
    min_s, max_e = 9e9, -1
    for clust in overlapping_call_clusts: 
        s,e = clust.get_med_start_end()
        min_s = min(min_s, s)
        max_e = max(max_e, e)
        all_starts_ends.append(s)
        all_starts_ends.append(e)
    
    all_starts_ends = sorted(np.unique(all_starts_ends))
    n_indivs = len(g.indivs)
    l = len(all_starts_ends)-1
    mus = np.zeros((n_indivs,l))
     
    for i in xrange(l):
        s, e = all_starts_ends[i], all_starts_ends[i+1] 
        X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        gX = g.GMM_genotype(X)
        print s, e, gX.n_clusts
        #Xs, s_idx_s, s_idx_e = g.get_sunk_gt_matrix(contig, s, e)
        #gXs = g.GMM_genotype(Xs)
        #g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, fn="./test/%d_%d.pdf"%(s, e))
        mus[:,i] = np.mean(X,1)    
    
    #print mus 
    #print np.corrcoef(mus)
    s, e = min_s, max_e
    X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
    gX = g.GMM_genotype(X)
    Xs, s_idx_s, s_idx_e = g.get_sunk_gt_matrix(contig, s, e)
    gXs = g.GMM_genotype(Xs)
    g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, fn="./test/%d_%d.png"%(s, e))
    c =  np.corrcoef(np.transpose(mus))
    print c
    return c


    
def get_best_gt(call, contig, g):
    
    max_bic_d = -9e9
    best_call = []
    best_se = []
    print call.starts
    print call.ends
    
    for s in np.unique(np.array(call.starts)):
        for e in np.unique(np.array(call.ends)):

            X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
            Xs, s_idx_s, s_idx_e = g.get_sunk_gt_matrix(contig, s, e)
            gX = g.GMM_genotype(X)
            gXs = g.GMM_genotype(Xs)

            d1 = gX.get_bic_delta()
            d2 = gXs.get_bic_delta()

            if max(d1, d2) > max_bic_d:
               max_bic_d = max(d1, d2) 
               best_inf = [gX, gXs, s,e,idx_s, idx_e, s_idx_s, s_idx_e]
            #if min(gX.min_bic, gXs.min_bic) < min_bic:
            #   min_bic = min(gX.min_bic, gXs.min_bic) 
            #   best_inf = [gX, gXs, s,e,idx_s, idx_e, s_idx_s, s_idx_e]
            g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, suffix="_TEST")
   
    gX, gXs, s, e, idx_s, idx_e, s_idx_s, s_idx_e = best_inf
    g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, suffix="_ZBEST")
    print s, e
    raw_input()
    

class filter_obj:
    def __init__(self, min_max_mu_d, max_overlap, max_mu_cp, singleton_min_sigma, dup_min_sigma, filter_X_linked=False, sex_lambda=None, p_thresh=1e-10):
        """
        min_max_mu_d: the minumum acceptible distance between gaussians for
        the maximum distance of any fit - ie, make sure that there is at least 
        one pair of adjacent guassians that looks at least this good

        max_mu_overlap: the maximum mean of the overlaps between
        adjacent fit gaussians to be accepted
    
        max_mu_cp: ignore any copy number variation in regions where the mean of all 
        the states is > max_mu_cp

        filter_X_linked: do a Fishers-exact test to see if there is male/female bias 
        in a call suggesting that X-chr paralogy is the culprit

        """
        self.min_max_mu_d = min_max_mu_d
        self.max_overlap = max_overlap
        self.max_mu_cp = max_mu_cp
        self.singleton_min_sigma = singleton_min_sigma
        self.dup_min_sigma = dup_min_sigma

        self.filter_X_linked = filter_X_linked
        self.p_thresh=p_thresh
        if sex_lambda:
            self.sex_lambda = sex_lambda
        else:
            self.sex_lambda = lambda x: x.rsplit("_",2)[-1]
    
    def is_X_linked(self, indivs, labels):
        u_labels = list(np.unique(labels))
        is_female = np.array([self.sex_lambda(indiv)=="F" for indiv in indivs])
        
        label_counts = [[label,np.sum(labels==label)] for label in u_labels]
        max_labs = [l[0] for l in sorted(label_counts, key=lambda x: -x[1])[:2]]
    
        #print "l\tfemale\tmale"
        rs=[]
        for label in max_labs:
            #print label, np.sum((labels==label)&is_female) , np.sum((labels==label)&(~is_female)) 
            rs.append([np.sum((labels==label)&is_female) , np.sum((labels==label)&(~is_female))])
        obs =  np.array(rs)

        #if np.amin(np.sum(obs,1))<3:
        #    return False

        odds, p = fisher_exact(obs)
        print "Fisher exact:", p

        if p<self.p_thresh:
            return True

        return False        
class GMM_gt(object):

    def __init__(self, X, gmm, labels, Z, params, bics, indivs):
        self.X = X
        self.mus = np.mean(self.X,1)
        self.gmm = gmm
        self.labels = labels
        self.Z = Z
        self.params = params
        self.bics = bics
        self.min_bic = min(self.bics)
        self.weights = self.gmm.weights
        
        self.indivs = indivs 
        self.n_clusts = np.unique(self.labels).shape[0]
        
        self.n_wnds = self.X.shape[1]

        self.f_correct = None 
        
        self.shaped_mus = np.reshape(self.mus, (self.mus.shape[0],1))
        self.l_probs, self.posterior_probs = gmm.eval(self.shaped_mus)
        
        #only in sklearn 0.14.1
        #self.score_samples = gmm.score_samples(X)


        #unique labels are the unique labels and uniq mus are the 
        #mean of the self.mus for each label
        
        self.label_to_mu = {}
        self.mu_to_labels = {}
        self.label_to_std = {}
        
        self.all_uniq_labels = []
        self.all_uniq_mus = []

        for l in np.unique(self.labels):
            self.label_to_mu[l] = np.mean(self.mus[self.labels==l])
            self.label_to_std[l] = np.std(self.mus[self.labels==l])
            self.mu_to_labels[np.mean(self.mus[self.labels==l])] = l
            self.all_uniq_labels.append(l)
            self.all_uniq_mus.append(np.mean(self.mus[self.labels==l]))

        self.all_uniq_labels = np.array(self.all_uniq_labels)
        self.all_uniq_mus = np.array(self.all_uniq_mus)

   
    def simple_plot(self, fn_out):

        cps = self.mus
        plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
        fig, ax_arr = plt.subplots(1,3)
        fig.set_figwidth(9)
        fig.set_figheight(4)
        axescolor  = '#f6f6f6'
        print ax_arr 
        print self.bics
        print self.params

        ax_arr[1].plot(self.params, self.bics)

        n, bins, patches = ax_arr[0].hist(cps,alpha=.9,ec='none',normed=1,color='#8DABFC',bins=len(cps)/10)
        #self.addGMM(gX.gmm, axarr[1,1], cps, gX.labels, overlaps)
        
        G_x=np.arange(0,max(cps)+1,.1)
        l = self.gmm.means.shape[0]
        
        for i in xrange(l):
            c = cm.hsv(float(i)/l,1)
            mu = self.gmm.means[i,0]
            var = self.gmm.covars[i][0][0]

            G_y = mlab.normpdf(G_x, mu, var**.5)*self.gmm.weights[i]
            ax_arr[0].plot(G_x,G_y,color=c)
            ax_arr[0].plot(mu,-.001,"^",ms=10,alpha=.7,color=c)
        
        if np.amax(cps)<2:
            ax_arr[0].set_xlim(0,2)
        ylims = ax_arr[0].get_ylim()
        if ylims[1] > 10:
            ax_arr[0].set_ylim(0,10)

        fig.sca(ax_arr[2]) 
        dendro = hclust.dendrogram(self.Z, orientation='right')
        ylims = ax_arr[2].get_ylim()
        ax_arr[2].set_ylim(ylims[0]-1, ylims[1]+1)

        fig.savefig(fn_out)

    def correct_order_proportion(self):
        #wnd_proportion_dir
        sorted_mus = np.sort(self.all_uniq_mus)
        labels = np.array(self.labels)
        if len(sorted_mus)==1: return 0.0
        
        t = 0 
        p = 0
        
        for i in xrange(sorted_mus.shape[0]-1):
            mu_0 = sorted_mus[i]
            l_0 = self.mu_to_labels[mu_0]
            
            mu_1 = sorted_mus[i+1]
            l_1 = self.mu_to_labels[mu_1]

            assert mu_0<mu_1
            
            wl_0 = np.where(labels==l_0)
            wl_1 = np.where(labels==l_1)
            t+=self.X.shape[1]
            """
            COULD traverse in...
            """
            #print self.X.shape
            #print self.X[wl_0,:]
            #print self.X[wl_1,:]
            #print wl_0[0].shape, wl_1[0].shape
            #print np.amax(self.X[wl_0,:], 1), np.amin(self.X[wl_1,:], 1)
            s=np.sum(np.amax(self.X[wl_0,:], 1)<np.amin(self.X[wl_1,:], 1))
            s2=np.sum(np.median(self.X[wl_0,:], 1)<np.median(self.X[wl_1,:], 1))
            p+=s
            #print "\t", mu_0, l_0, mu_1, l_1, s, self.X.shape[1], s/float(self.X.shape[1]),s2/float(self.X.shape[1])      
        
        return float(p)/t
    
    def get_min_z_dist(self):

        min_z_dist = 9e9

        sorted_mus = np.sort(self.all_uniq_mus)
        labels = np.array(self.labels)
        
        delta = np.diff(sorted_mus)
        std_lefts = []
        std_rights = []

        for i in xrange(sorted_mus.shape[0]-1):
            mu_left = sorted_mus[i]
            mu_right = sorted_mus[i+1]
            
            l_left = self.mu_to_labels[mu_left]
            l_right = self.mu_to_labels[mu_right]
            
            std_left = max(1e-9, self.label_to_std[l_left])
            std_right = max(1e-9, self.label_to_std[l_right])
            
            std_lefts.append(std_left)
            std_rights.append(std_right)
        
        std_lefts = np.array(std_lefts)
        std_rights = np.array(std_rights)

        min_l = np.amin(np.absolute(delta/std_lefts)) 
        min_r = np.amin(np.absolute(delta/std_rights)) 
        
        return min(min_l, min_r)


    def fail_filter(self, filt):
        
        if self.n_clusts ==1 or self.gmm.means.shape[0] == 1:  
            return True

        u_o, med_o, overlaps = assess_GT_overlaps(self.gmm)
        max_overlap_stat = sorted(overlaps, key = lambda x: max(x['os']))[-1]
        
        """
        filter by the fractional overlap between adjacent clusters
        """
        if max(max_overlap_stat['os'])>=filt.max_overlap:
            return True

        mu_mu_d, min_mu_d, max_mu_d = self.get_mean_min_max_inter_mu_dist()
        
        """ 
        filter by the maximum distance between clusters, 
        it must be greater than min
        """
        if (max_mu_d < filt.min_max_mu_d):
            return True

        """
        filter out regions with really high copy
        """
        if np.mean(self.mus)>filt.max_mu_cp:
            return True
        
        """
        if flagged, then test for x-linked association and remove  
        """
        if filt.filter_X_linked and filt.is_X_linked(self.indivs, self.labels):
            return True
        
        """
        if singleton, impose strict theshold
        """
        if  ( self.n_clusts == 2 ) and ( min(np.sum(self.labels==l) for l in np.unique(self.labels)) == 1 ):
            min_z = self.get_min_z_dist()
            if min_z < filt.singleton_min_sigma:
                return True
        
        """
        if a dup, force stricter genotype filtering
        gts_by_indiv, gts_to_label, labels_to_gt = gX.get_gts_by_indiv()
        """
        if (np.amax(self.mus)>2.5):
            min_z = self.get_min_z_dist()
            if min_z < filt.dup_min_sigma:
                return True
            
        return False
        
        #mean_mu_delta = self.get_mean_inter_mu_dist()
        
        #if self.f_correct == None: 
        #    self.f_correct = self.correct_order_proportion()
        #if (self.f_correct <= frac_dir_min):
        #    return True


    def output_filter_data(self, F_filt, contig, s, e):
        
        if self.f_correct == None:
            self.f_correct = self.correct_order_proportion()
        
        mu_mu_d, min_mu_d, max_mu_d = self.get_mean_min_max_inter_mu_dist()

        n_clusts = self.n_clusts 
        min_AC = min(np.sum(self.labels==l) for l in np.unique(self.labels))

        n_wnds = self.n_wnds
        min_z = self.get_min_z_dist()
        
        bic_delta = self.get_bic_delta() 
        
        F_filt.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(contig, s, e, 
                                                                         mu_mu_d, 
                                                                         max_mu_d, 
                                                                         min_mu_d, 
                                                                         self.f_correct,
                                                                         min_z, 
                                                                         n_wnds,
                                                                         bic_delta,
                                                                         n_clusts,
                                                                         min_AC))
    
    def get_mean_min_max_inter_mu_dist(self):
        s_mus = np.sort(self.all_uniq_mus)
        ds = np.diff(s_mus)
        if ds.shape[0] == 0:
            return 0, 0, 0

        return np.mean(ds), np.amin(ds), np.amax(ds) 
        
    def get_gts_by_indiv(self):
        
        cp_2_thresh=1.0
        mode_label = int(mode(self.labels)[0][0])
        
        indiv_labels = np.array(self.labels)
        
        mu_args = np.argsort(self.all_uniq_mus) 
          
        ordered_labels = self.all_uniq_labels[mu_args]
        ordered_mus = self.all_uniq_mus[mu_args] 
        d_from_2 = np.absolute(ordered_mus-2.0)
        
        labels_to_gt = {}
        """
        if there is something that looks like a 2, use it to callibrate others
        assign 2 to the closest 1, then assign the rest as +-1 in the order
        make sure that you aren't assigning -1 genotypes
        then finally, consolidate w/ the mus 
        """
        if np.amin(d_from_2)<cp_2_thresh:
            idx = np.argmin(d_from_2)
            idx_cp = 2
        else:
            idx = 0
            idx_cp = round(ordered_mus[0])

        for i,l in enumerate(ordered_labels): 
            labels_to_gt[l] = idx_cp-(idx-i)
        
        ## ensure no -1s
        while min(labels_to_gt.values())<0:
            print "<0's detected..."
            new_labels_to_gt = {}
            for l, gt in labels_to_gt.iteritems():
                new_labels_to_gt[l] = gt+1
            labels_to_gt = new_labels_to_gt
       
        ##correct for odd major alleles out of HWE 
        if (labels_to_gt[mode_label] %2 == 1) and np.sum(indiv_labels==mode_label) >= .5*(indiv_labels.shape[0]):
            d=0
            if self.label_to_mu[mode_label]-labels_to_gt[mode_label]>0 or min(labels_to_gt.values())==0:
                d=1
            else:
                d=-1
            new_labels_to_gt = {}
            for l, gt in labels_to_gt.iteritems():
                new_labels_to_gt[l] = gt+d
            labels_to_gt = new_labels_to_gt
        
        gts_by_indiv = {}
        for i, indiv in enumerate(self.indivs):  
            gts_by_indiv[indiv] = int(labels_to_gt[self.labels[i]]) 
        
        new_labels_to_gt = {k:int(v) for k,v in labels_to_gt.iteritems()}
        labels_to_gt = new_labels_to_gt

        gt_to_labels = {v:k for k,v in labels_to_gt.iteritems()} 

        return gts_by_indiv, gt_to_labels, labels_to_gt
            
    def get_gt_lls_by_indiv(self, labels_to_gt):
        

        gt_lls_by_indiv = {}

        for i, indiv in enumerate(self.indivs):  
            lls = {gt:max(np.log(self.posterior_probs[i,label])/np.log(10),-1000) for label, gt in labels_to_gt.iteritems() }
            gt_lls_by_indiv[indiv] = lls 

        return gt_lls_by_indiv


    def get_bic_delta(self):
        
        if self.n_clusts+1 in self.params:
            idx_bic_r = self.params.index(self.n_clusts+1) 
            bic_r = self.bics[idx_bic_r]
        else:
            bic_r = self.min_bic
    
        if self.n_clusts-1 in self.params:
            idx_bic_l = self.params.index(self.n_clusts-1) 
            bic_l = self.bics[idx_bic_l]
        else:
            bic_l = self.min_bic
        
        delta = bic_r-self.min_bic+bic_l-self.min_bic
        return delta

    def get_cp(self, indiv, g):
        
        idx = g.indivs.index(indiv_id)
        

    def is_var(self, indiv_id, g, force_not_mode = False):

        idx = g.indivs.index(indiv_id)
        
        if (not force_not_mode) and len(np.unique(self.labels))>2:
            return True
        else:
            idx = g.indivs.index(indiv_id)
            m = mode(self.labels)[0]
            if self.labels[idx] != m:
                return True
        
        return False

def eval_G(G, x):

    u, v = G
    s = np.sqrt(v)
    sq2pi = np.power(2*np.pi,0.5)

    y = (1/(s*sq2pi)) * np.exp( -1*((x-u)*(x-u))/(2*s*s) )
    return y

def get_intersection(G1, G2, ws, tol=0.01):
    """
    Get the integral of the shared area between two gaussians
    the intersection is computed exactly, but the integral is computed numerically
    """
    tol = min(ws)/1000
    #sort so G1.mu < G2.mu
    #ui < uj
    oGs = [G1, G2] 
    ows = ws
    Gs, ws = [], []
    args = np.argsort([G1[0],G2[0]])
    
    for i in args:
        Gs.append(oGs[i])
        ws.append(ows[i])
    ui, vi = Gs[0]
    uj, vj = Gs[1]
    si, sj = np.sqrt(vi), np.sqrt(vj)
    al, be = ws
    
    if si == sj:
        x=(ui+uj)/2.0
    else:
        sq2pi = np.power(2*np.pi,0.5)
        c = (2*si*si*sj*sj) * ( np.log( al/(si*sq2pi) ) - np.log( be/(sj*sq2pi) ) )
        c = c  + (si*si*uj*uj)-(sj*sj*ui*ui)
        b = -((2*uj*si*si)-(2*ui*sj*sj))
        a = (si*si)-(sj*sj)
        
        q=(b**2 - 4*a*c)
        if q<0: 
            x=None
        else:
            x1 = (-b + np.sqrt(q)) / (2*a)
            x2 = (-b - np.sqrt(q)) / (2*a)
            
            x=x1
            if (x1 < ui and x1 < uj) or (x1 > ui and x1 > uj):
                x=x2
    
    if x==None:
        return None, None, 1, 1 

    y = al*eval_G(G1, x) 
    
    k=5
    mn = ui - k*si
    mx = uj + k*sj

    o_gran_mx = (mx-x)/tol
    o_gran_mn = (x-mn)/tol
    
    while np.absolute((mx-x)/tol)>1e7 and k>2:
        mn = ui - k*si
        mx = uj + k*sj
        k-=.2
    
    if o_gran_mx!=(mx-x)/tol:
        print "\tgaussian intercept granularity reduction:", o_gran_mx, (mx-x)/tol, k
        print "\tgaussian intercept granularity reduction:", o_gran_mn, (x-mn)/tol, k

    if np.absolute((mx-x)/tol)>1e7:
        return None, None, 1, 1 


    xis = np.arange(x,mx, tol)
    xjs = np.arange(mn,x, tol)

    i_integral = np.sum(mlab.normpdf(xis, ui, si)*al)*tol
    j_integral = np.sum(mlab.normpdf(xjs, uj, sj)*be)*tol
    overlap = i_integral+j_integral

    return x, y, overlap/al, overlap/be

def assess_GT_overlaps(gmm):

    overlaps = []

    l = gmm.means.shape[0] 
    us = []
    ss = []
    ws = []
    for i in xrange(l):
        u = gmm.means[i,0]
        s = gmm.covars[i][0][0]**.5
        w = gmm.weights[i]
        if w==0: continue
        us.append(u)
        ss.append(s)
        ws.append(w)
    
    sort_mu_args = np.argsort(np.array(us))
    all_os = []
    for k in xrange(len(sort_mu_args)-1):
        i, j = sort_mu_args[k], sort_mu_args[k+1] 
        u_1, u_2 = us[i], us[j]
        s1, s2 = ss[i], ss[j]
        G1, G2 = [u_1,s1*s1], [u_2,s2*s2]
        w1, w2 = ws[i], ws[j]
        t = w1+w2
        w1, w2 = w1/t, w2/t
        x, y, o1, o2 = get_intersection(G1, G2, [w1,w2], tol=0.01)
        s_dist=[np.absolute(u_1-u_2)/s1,np.absolute(u_1-u_2)/s2]
        all_os+=[o1,o2]
        overlaps.append({"us":tuple([u_1,u_2]),"os":tuple([o1, o2]),"ss":tuple([s1,s2]), "ws":tuple([w1,w2]), "sdist":tuple(s_dist)})
    
    u_o, med_o = np.mean(all_os), np.median(all_os)
    return u_o, med_o, overlaps
    
        
def output(g, contig, s, e, F_gt, F_call, F_filt, F_VCF, filt, include_indivs=None, plot=False, v=False):

    print "%s %d %d"%(contig, s, e)
    stdout.flush()

    if include_indivs!=None and len(include_indivs) == 1:
        return
    
    X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
    gX = g.GMM_genotype(X, include_indivs = include_indivs)
    u_o, med_o, overlaps = assess_GT_overlaps(gX.gmm)

    if gX.fail_filter(filt):
        print "***********FAILED************"
    if gX.n_clusts ==1:  
        print "***********1_CLUST************"

    if gX.n_clusts == 1 or gX.fail_filter(filt):
        return

    F_call.write("%s\t%d\t%d\n"%(contig, s, e))
    g.output(F_gt, F_VCF, gX, contig, s, e, v=v)

    gX.output_filter_data(F_filt, contig, s, e)
    
     
    if plot:
        print "plotting %s %d %d"%(contig, s, e)
        Xs, s_idx_s, s_idx_e = g.get_sunk_gt_matrix(contig, s, e)
        gXs = g.GMM_genotype(Xs, include_indivs = include_indivs)
        g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, overlaps, fn="./plotting/test/%s_%d_%d.png"%(contig, s, e))


class genotyper:
    
    def setup_output(self, FOUT, FFILT, F_VCF):
        outstr = "contig\tstart\tend\t%s\n"%("\t".join(self.indivs))
        FOUT.write(outstr)
        
        VCF_outstr = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"%("\t".join(self.indivs))
        F_VCF.write(VCF_outstr)
        
        FFILT.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format("chr", 
                                                                 "start", 
                                                                 "end", 
                                                                 "mu_mu_d", 
                                                                 "max_mu_d", 
                                                                 "min_mu_d", 
                                                                 "f_correct_direction",
                                                                 "min_z",
                                                                 "wnd_size", 
                                                                 "bic_delta",
                                                                 "n_clusts",
                                                                 "min_allele_count"))

    def output(self, FOUT, F_VCF, gX, contig, s, e, v=False):
        
        outstr = "%s\t%d\t%d"%(contig, s, e)
        gts_by_indiv, gts_to_label, labels_to_gt = gX.get_gts_by_indiv()
        
        cns = sorted(gts_to_label.keys())
        cn_range=range(0,max(cns)+1)
        
        gt_lls_by_indiv = gX.get_gt_lls_by_indiv(labels_to_gt)
        
        
        #reference hap_cn is already 1
        hap_cns = [1]
        if 0 in cns or 1 in cns:
            hap_cns = [1,0]
        
        if max(cns)>2:
            max_hap = int(math.ceil(max(cns)/2.0))
            hap_cns += range(2,max_hap+1) 
        
        hap_cns_to_labels = {1:0}
        
        i=1
        for hcn in hap_cns:
            if hcn==1: continue
            hap_cns_to_labels[hcn]=i
            i+=1

        #cns = 0,1,2 
        #F(j/k) = (k*(k+1)/2)+j 
        GTs_to_cps = {} 
        cps_to_GTs = {} 
        jk_tups = []
        jk_tups_to_cp = {} 
        for i in xrange(len(hap_cns)):
            for j in xrange(i,len(hap_cns)):
                c1, c2 = hap_cns[i], hap_cns[j]
                cp = c1+c2
                lbl = "%d/%d"%(hap_cns_to_labels[c1], hap_cns_to_labels[c2])
                GTs_to_cps[lbl] = cp 
                jk_tup = tuple([hap_cns_to_labels[c1], hap_cns_to_labels[c2]])
                jk_tups.append(jk_tup)
                jk_tups_to_cp[jk_tup] = cp

                if not cp in cps_to_GTs: 
                    cps_to_GTs[cp] = []

                cps_to_GTs[cp].append(lbl)
        
        if v:
            print "cps_to_GTs", cps_to_GTs
            print "jk_tups_to_cp", jk_tups_to_cp
            print "GTs_to_cps:", GTs_to_cps

        """
        VCF OUTPUT
        GT:CN:CNL:GL:PL
        CN = copy number
        CNP = Copy Number likelihood

        #PL.. phred scaled genotype likelihood...? MAYBE?
        """
        V_outstr = "{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}"
        ALTS = ",".join(["<CN%d>"%c for c in hap_cns if c!=1])

        VCF_contig = contig.replace("chr","")

        if ALTS=="<CN2>":
            SVTYPE="DUP"
        elif ALTS=="<CN0>":
            SVTYPE="DEL"
        else:
            SVTYPE="CNV"
        
        INFO="END=%d;SVTYPE=%s;LPROBS=%.2f"%(e, SVTYPE ,np.sum(gX.l_probs))
        V_outstr = V_outstr.format(CHROM=VCF_contig,
                                   POS=s+1,
                                   ID="%s_%d_%d"%(contig, s+1, e),
                                   REF=self.fasta['%s:%d'%(VCF_contig,s+1)],
                                   ALT=ALTS,
                                   QUAL='.',
                                   INFO=INFO,
                                   FILTER="PASS",
                                   FORMAT="GT:CN:GL:PL:CNL")
        
        V_data = []
        ordered_cps = []
        for indiv in self.indivs:
            if indiv in gts_by_indiv:
                ordered_cps.append(gts_by_indiv[indiv])
                
                s_CNL=",".join(["%.2f"%(cn not in cns and -1000 or gt_lls_by_indiv[indiv][cn]) for cn in cn_range])
                s_GLs = ""
                COPY=gts_by_indiv[indiv]
                GT = cps_to_GTs[COPY][0]
                
                GLs = np.zeros(len(jk_tups))
                PLs = np.zeros(len(jk_tups))
                #F(j/k) = (k*(k+1)/2)+j 
                for jk in jk_tups:
                    j,k = jk
                    p = (k*(k+1))/2+j
                    cp = jk_tups_to_cp[jk]
                    if not cp in cns:
                        GLs[p] = -1000
                        PLs[p] = 1000
                    else:
                        GLs[p] = gt_lls_by_indiv[indiv][cp] - (np.log(len(cps_to_GTs[cp]))  / np.log(10))
                        PLs[p] = GLs[p] == -1000 and 1000 or int(round(-10*GLs[p]))
                        
                sGLs = ",".join(["%.2f"%l for l in GLs])
                sPLs = ",".join(["%d"%l for l in PLs])
                s="{GT}:{COPY}:{GLs}:{PLs}:{CNL}".format(GT=GT,
                                                         COPY=COPY,
                                                         GLs = sGLs,
                                                         PLs=sPLs,
                                                         CNL=s_CNL)
                V_data.append(s)
            else:
                ordered_cps.append(-1)
                V_data.append("./.:0:0:0")

        outstr = "%s\t%s\n"%(outstr, "\t".join("%d"%gt for gt in ordered_cps))
        V_outstr  = "%s\t%s\n"%(V_outstr, "\t".join(V_data))
        
        FOUT.write(outstr)
        F_VCF.write(V_outstr)
        
        if v:
            print V_outstr
            print gX.l_probs
            print gX.posterior_probs
            print gX.labels
            print np.unique(gX.labels)
            print "gts_to_label:", gts_to_label, "labels_to_gt:", labels_to_gt
             
    def init_on_indiv_DTS_files(self, **kwargs):

        g = gglob.init_from_DTS(**kwags)

        self.indivs = g.indivs
        self.wnd_starts = g.wnd_starts
        self.wnd_ends = g.wnd_ends
        self.sunk_wnd_starts = g.sunk_wnd_starts
        self.sunk_wnd_ends = g.sunk_wnd_ends
        
        self.cp_matrix = g.cp_matrix
        self.sunk_cp_matrix = g.sunk_cp_matrix

    
    def init_on_gglob(self, contig, fn_gglob, subset_indivs):
        
        g = gglob.init_from_gglob_dir(contig, fn_gglob, indiv_subset=subset_indivs)
        
        self.indivs = g.indivs
        self.wnd_starts = g.wnd_starts
        self.wnd_ends = g.wnd_ends
        self.sunk_wnd_starts = g.sunk_wnd_starts
        self.sunk_wnd_ends = g.sunk_wnd_ends
        
        self.cp_matrix = g.cp_matrix
        self.sunk_cp_matrix = g.sunk_cp_matrix
        

    def __init__(self, contig, **kwargs): 

        self.gglob_dir = kwargs.get("gglob_dir", None) 
        self.plot_dir  = kwargs.get("plot_dir", None)
        
        subset_indivs  = kwargs.get("subset_indivs", None)
        fn_fasta  = kwargs.get("fn_fa", None)
        

        self.contig = contig
        self.indivs = []
        self.wnd_starts = None
        self.wnd_ends = None
        self.sunk_wnd_starts = None
        self.sunk_wnd_ends = None
        
        self.cp_matrix = None
        self.sunk_cp_matrix = None
        
        if self.gglob_dir:
            self.init_on_gglob(self.gglob_dir, self.contig, subset_indivs) 
        else:
            self.init_on_indiv_DTS_files(self, **kwargs)
    
        k = self.cp_matrix.shape[0]
        print >>stderr, "loading %d genomes..."%(k)
        t = time.time()
        print >>stderr, "done (%fs)"%(time.time()-t)
        
        if fn_fasta:
            self.fasta = FastaHack(fn_fasta)
       
    def addGMM(self, gmm, ax, X, labels, overlaps=None):
        
        G_x=np.arange(0,max(X)+1,.1)
        l = gmm.means.shape[0]
        u_labels = np.unique(labels)
        print u_labels

        for i in xrange(l):
            #if not i in u_labels: continue
            c = cm.hsv(float(i)/l,1)
            mu = gmm.means[i,0]
            var = gmm.covars[i][0][0]

            print mu, var, var**.5
            
            G_y = mlab.normpdf(G_x, mu, var**.5)*gmm.weights[i]
            ax.plot(G_x,G_y,color=c)
            ax.plot(mu,-.001,"^",ms=10,alpha=.7,color=c)
        
        if overlaps!=None:
            ymax=ax.get_ylim()[1]
            y=ymax-.2
            xmin=ax.get_xlim()[0]

            all_os = []
            for o in overlaps:
                y-=.4
                us = o['us']
                os = o['os']
                ss = o['ss']
                ws = o['ws']

                x=(us[0]+us[1])/2.0
                d=us[0]-us[1]
                ds1 = d/ss[0]
                ds2 = d/ss[1]

                o1, o2 = os[0], os[1]
                if o1 == None: o1, o2 = 1.0, 1.0
                all_os+=[o1, o2]
                
                ax.text(x,y,"%.2f %.2f %.2f %.2f %.2f"%(o1, ds1, us[0], ss[0], ws[0] ), fontsize=6, horizontalalignment='center', verticalalignment='center')
                ax.text(x,y-.15,"%.2f %.2f %.2f %.2f %.2f"%(o2, ds2, us[1], ss[1], ws[1] ), fontsize=6, verticalalignment='center', horizontalalignment='center')
             
            ax.text(xmin+1,ymax-.2,"%.2f %.2f"%(np.mean(all_os), np.max(all_os)), fontsize=8, verticalalignment='center', horizontalalignment='right')
            
    def aug_dendrogram(self, ax, ddata):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            ax.plot(x, y, 'ro')
            ax.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    def plot(self, gX, gXs, chr, start, end, idx_s, idx_e, s_idx_s, s_idx_e, overlaps, fn="test_gt.pdf"):
        
        X = gX.X
        Xs = gXs.X

        cps = np.mean(X, 1)
        sunk_cps = np.mean(Xs, 1)
        
        plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
        fig, axarr = plt.subplots(3, 3)
        fig.set_figwidth(11)
        fig.set_figheight(8.5)
        axescolor  = '#f6f6f6'
        
        axarr[0,0].plot(cps, sunk_cps, 'ro', alpha=0.2)
        axarr[0,0].set_xlim(-0.10,max(cps)+1)
        axarr[0,0].set_ylim(-0.10,max(sunk_cps)+1)
        

        if len(cps)<=50:
            b = len(cps)/2
        else:
            b = len(cps)/10

        n, bins, patches = axarr[1,1].hist(cps,alpha=.9,ec='none',normed=1,color='#8DABFC',bins=b)
        self.addGMM(gX.gmm, axarr[1,1], cps, gX.labels, overlaps)
        
        fig.sca(axarr[0,2]) 
        dendro = hclust.dendrogram(gX.Z, orientation='right')
        ylims = axarr[0,2].get_ylim()
        axarr[0,2].set_ylim(ylims[0]-1, ylims[1]+1)
        #self.aug_dendrogram(axarr[0,2], dendro)
        
        if np.sum(np.isnan(sunk_cps)) == 0:
            n, bins, patches = axarr[1,0].hist(sunk_cps,alpha=.9,ec='none',normed=1,color='#FCDE8D',bins=b)
            self.addGMM(gXs.gmm, axarr[1,0], sunk_cps, gXs.labels)
            axarr[0,1].plot(gX.params, gX.bics, 'ro', ms=5)
            axarr[0,1].plot(gXs.params, gXs.bics, 'go', ms=5)
            amin = np.argmin(np.array(gX.bics))
            axarr[0,1].plot(gX.params[amin], gX.bics[amin], 'bo', ms=7)
        
        fig.sca(axarr[1,2]) 
        dendro = hclust.dendrogram(gXs.Z, orientation='right')
        ylims = axarr[1,2].get_ylim()
        axarr[1,2].set_ylim(ylims[0]-1000, ylims[1]+1000)
        #self.aug_dendrogram(axarr[1,2], dendro)

        #plot actual position

        #def get_gt_matrix(self, contig, start, end, vb=False):
        #    assert contig == self.contig
        #    start_idx = np.searchsorted(self.wnd_starts, start)
        #    end_idx = np.searchsorted(self.wnd_ends, end)
        
        #idx_s, idx_e = np.where(self.wnd_starts==start)[0], np.where(self.wnd_ends==end)[0]+1
        #if idx_e <= idx_s: idx_e = idx_s+1
        #s_idx_s, s_idx_e = np.searchsorted(self.sunk_wnd_starts, start),  np.searchsorted(self.sunk_wnd_ends, end)+1
        #if s_idx_e <= s_idx_s: s_idx_e = s_idx_s+1
        
        xs = (self.wnd_starts[idx_s:idx_e]+self.wnd_ends[idx_s:idx_e])/2.0
        s_xs = (self.sunk_wnd_starts[s_idx_s:s_idx_e]+self.sunk_wnd_ends[s_idx_s:s_idx_e])/2.0
        #print "shapes", X.shape, Xs.shape

        for i in xrange(X.shape[0]):
            axarr[2,1].plot(xs, X[i,:])
            axarr[2,0].plot(s_xs, Xs[i,:])
        axarr[2,1].set_xlim(start,end) 
        axarr[2,0].set_xlim(start,end) 
        #fig.savefig("%s/%s-%d-%d%s.png"%(self.plot_dir, chr, start, end, suffix))
        fig.savefig(fn)
        plt.close()

        
    def get_gt_matrix(self, contig, start, end, vb=False):
        assert contig == self.contig
        
        start_idx = np.searchsorted(self.wnd_starts, start)
        end_idx = np.searchsorted(self.wnd_ends, end)
        
        if end_idx<=start_idx:
            end_idx = start_idx+1
        elif end_idx-start_idx>1:
            start_idx+=1
        
        X = self.cp_matrix[:,start_idx:end_idx]
         
        return X, start_idx, end_idx
    
    def get_sunk_gt_matrix(self, contig, start, end, vb=False):
        assert contig == self.contig
        
        start_idx = np.searchsorted(self.sunk_wnd_starts, start)
        end_idx = np.searchsorted(self.sunk_wnd_ends, end)
        
        #print "BEGIN:", start_idx, end_idx 
        if end_idx<=start_idx:
            end_idx = start_idx+1
        elif end_idx-start_idx>1:
            start_idx+=1
        #print "END:", start_idx, end_idx 

        X = self.sunk_cp_matrix[:,start_idx:end_idx]
         
        return X, start_idx, end_idx
    
    def get_gt_matrix_mu(self, contig, start, end):
        assert contig == self.contig
        start_idx = np.searchsorted(self.wnd_starts, start)
        end_idx = np.searchsorted(self.wnd_ends, end)+1
        X = self.cp_matrix[:,start_idx:end_idx]
        
        return X

    def s_score(self, X, labels):
        return metrics.silhouette_score(X, labels) 
    
    def fit_GMM(self, X, init_means, init_vars, init_weights, n_iter=1000):
        #min_covar=1e-5
        min_covar=1e-4
        n_components = len(init_means)
        #gmm = mixture.GMM(n_components, 'spherical')
        #gmm = mixture.GMM(n_components, 'diag')
        gmm = mixture.GMM(n_components, 'spherical', min_covar=min_covar)
        gmm.means = np.reshape(np.array(init_means),(len(init_means),1))
        gmm.weights = np.array(init_weights)
        
        #covars = np.array([v[0][0] for v in gmm.covars])
        #covars = np.array([np.reshape(np.array([max(v, min_covar)]),(1,1)) for v in init_vars])
        covars = np.array([max(v, min_covar)for v in init_vars])
        gmm.covars = covars
        
        #gmm.fit(X, n_iter=n_iter, init_params='c')
        #gmm.fit(X, n_iter=n_iter, init_params='cmw')
        #gmm.fit(X, n_iter=n_iter, init_params='c')
        gmm.fit(X, n_iter=n_iter, init_params='')
        
        labels = gmm.predict(X)
        
        bic = -2*gmm.score(X).sum() + (3*n_components)*np.log(X.shape[0])
        aic = -2*gmm.score(X).sum() + 2*(3*n_components)
        
        #print np.unique(labels)
        #for i,l in enumerate(list(np.unique(labels))):
        #    print init_means[i], init_vars[i], init_weights[i], gmm.score(X)[labels==l].sum() 
        #print bic

        #pdb.set_trace()
        return gmm, labels, aic 
    

    def GMM_genotype(self, X, include_indivs = None, FOUT = None, overload_indivs = None):
        """
        GMM genotyping
        merge_overlap_thresh, if -1, don't ever merge, however, 
        otherwise, if overlap > merge_overlap_thresh, then merge and recalculate 
        
        include indivs: only genotype a subset of indivs in the gt object

        overlaod indivs: this is if the gt object has nothing to do with the X you are
                         passing to the function
        """
        if include_indivs:
            l = len(include_indivs)
            new_X = np.zeros((l, X.shape[1]))
            
            for i, indiv in enumerate(include_indivs):
                j = self.indivs.index(indiv)
                new_X[i] = X[j] 
                
            X=new_X
       
        mus = np.mean(X,1)
        mus = np.reshape(mus, (mus.shape[0],1))
        dist_mat = dist.pdist(mus)
        print "h_clustering..."
        t = time.time()
        Z = hclust.linkage(mus, method='centroid', metric='euclidean')
        print "done %fs"%(time.time()-t)
        params, bics, gmms, all_labels = [], [], [], []
        
        print "assessing genotypes" 
        t = time.time()

        prev_grps = np.array([])
        for k in np.arange(.2, 1.5,  0.001):
            grps = hclust.fcluster(Z, k, criterion='distance')
            if np.all(grps == prev_grps): continue
            
            init_mus, init_vars, init_weights = self.initialize(mus, grps) 
            
            if len(init_mus)>30 and len(bics)>0: continue
            #pdb.set_trace()
            
            gmm, labels, ic = self.fit_GMM(mus, init_mus, init_vars, init_weights)

            params.append(len(init_mus))
            bics.append(ic)
            gmms.append(gmm)
            all_labels.append(labels)
            prev_grps = grps 
            
        print "done %fs"%(time.time()-t)
        idx = np.argmin(bics)
        #print params
        #print np.where(np.array(params)==3)
        #idx = np.where(np.array(params)==3)[0][0]
        gmm = gmms[idx]
        labels = all_labels[idx]
        

        ####NOW, finally merge calls that are too close 
        n_labels = np.unique(labels).shape[0] 
        if n_labels>1:
            gmm, labels = self.final_call_merge(gmm, labels, mus) 

        if include_indivs == None: 
            include_indivs = self.indivs
        
        if overload_indivs != None:
            return GMM_gt(X, gmm, labels, Z, params, bics, overload_indivs)
        else:
            return GMM_gt(X, gmm, labels, Z, params, bics, include_indivs)

    def final_call_merge(self, gmm, labels,mus, max_overlap=0.5, min_dist=0.6):
        """
        take the final min_bic call and merge calls that are too close   
        """
        #max_overlap=8

        u_o, med_o, overlaps = assess_GT_overlaps(gmm)
        max_overlap_stat = sorted(overlaps, key = lambda x: max(x['os']))[-1]
        n_labels = np.unique(labels).shape[0] 
        #pdb.set_trace()
        while (max(max_overlap_stat['os']) > max_overlap) and n_labels>1:
            u1, u2 = max_overlap_stat['us'] 
            l1, l2 = np.where(gmm.means==u1)[0][0], np.where(gmm.means==u2)[0][0]
            labels[labels==l2] = l1
            
            init_mus, init_vars, init_weights = self.initialize(mus, labels) 
            gmm, labels, ic = self.fit_GMM(mus, init_mus, init_vars, init_weights, n_iter=1000)

            n_labels = np.unique(labels).shape[0] 
            if n_labels>1:
                u_o, med_o, overlaps = assess_GT_overlaps(gmm)
                max_overlap_stat = sorted(overlaps, key = lambda x: max(x['os']))[-1]

        n_labels = np.unique(labels).shape[0] 
        #pdb.set_trace()

        if n_labels==1:
            return gmm, labels
        
        u_o, med_o, overlaps = assess_GT_overlaps(gmm)
        max_overlap_stat = sorted(overlaps, key = lambda x: np.absolute(x['us'][0]-x['us'][1]))[0]
        d = np.absolute(max_overlap_stat['us'][0]- max_overlap_stat['us'][1])
        #note, this is the min_sd of the max_overlap... 
        #maybe this is always the min_sd? but, perhaps not??
        sd = min(max_overlap_stat['sdist'])

        
        print overlaps
        while (d < min_dist ) and n_labels>1:
            u1, u2 = max_overlap_stat['us'] 
            l1, l2 = np.where(gmm.means==u1)[0][0], np.where(gmm.means==u2)[0][0]
            labels[labels==l2] = l1
            
            init_mus, init_vars, init_weights = self.initialize(mus, labels) 
            gmm, labels, ic = self.fit_GMM(mus, init_mus, init_vars, init_weights, n_iter=1000)

            u_o, med_o, overlaps = assess_GT_overlaps(gmm)
            print overlaps
            n_labels = np.unique(labels).shape[0] 

            if n_labels>1:
                max_overlap_stat = sorted(overlaps, key = lambda x: np.absolute(x['us'][0]-x['us'][1]))[0]
                d = np.absolute(max_overlap_stat['us'][0]- max_overlap_stat['us'][1])
                sd = min(max_overlap_stat['sdist'])
        
        return gmm, labels

    def initialize(self, X, grps):

        uniqs = np.unique(grps)
        init_mus =  []
        init_weights =  []
        init_vars =  []
        
        l = X.shape[0] 
        for grp in uniqs:
            init_mus.append(np.mean(X[grps==grp]))
            init_vars.append(np.var(X[grps==grp]))
            init_weights.append(float(np.sum(grps==grp))/l)
        
        return init_mus, init_vars, init_weights

    def ms_genotype(self, X):
        """
        mean shift based genotyping
        nlogn -> n^2
        """
        bandwidth = cluster.estimate_bandwidth(X, quantile=0.3)
        ms = cluster.MeanShift(bandwidth=bandwidth, bin_seeding=True)
        try:
            ms.fit(X)
        except:
            return -1, -1, -1

        labels = ms.labels_
        n_classes = np.shape(np.unique(labels))[0]
        if n_classes == 1: 
            s=-1
        else:
            s = self.s_score(X, labels)
        
        return labels, n_classes, s
    
    def AP_genotype(self, X): 
        """
        affinity propegation based genotyping
        n^2 in points
        """
        #ap = cluster.AffinityPropagation(damping=.9, preference=-200)
        ap = cluster.AffinityPropagation(damping=0.99)
        S = metrics.euclidean_distances(X)
        ap.fit(1-S)
        
        labels =  ap.labels_
        n_classes = np.shape(np.unique(labels))[0]
        if n_classes == 1: 
            s=-1
        else:
            s = self.s_score(X, labels)
        
        return labels, n_classes, s
    
    def DBS_genotype(self, X):
        """
        DBScan based genotyping
        likely bad...
        """
        dbs = cluster.DBSCAN(eps=1, min_samples=3)    
        dbs.fit(X)
        
        labels =  dbs.labels_
        n_classes = np.shape(np.unique(labels))[0]
        if n_classes == 1: 
            s=-1
        else:
            s = self.s_score(X, labels)
        
        return labels, n_classes, s
    
