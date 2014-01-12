import glob
from sys import stderr
import numpy as np


import time
import matplotlib.pyplot as plt
import matplotlib.colors as mCols
import matplotlib.cm as cm
import matplotlib.mlab as mlab

##local repo
from sklearn import cluster 
from sklearn import metrics
from sklearn import mixture


from sets import Set

import math
import random
from scipy.stats import norm

import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hclust

from wnd_cp_data import wnd_cp_indiv
from gglob import gglob

import cluster as m_cluster

from scipy.stats.mstats import mode

import IPython



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




def get_correlation_matrix(starts_ends, g, contig):

    n_indivs = len(g.indivs)
    l = len(starts_ends)-1
    mus = np.zeros((n_indivs,l))

    starts_ends = sorted(np.unique(starts_ends))

    n_indivs = len(g.indivs)
    l = len(starts_ends)-1
    mus = np.zeros((l,n_indivs))
     
    for i in xrange(l):
        s, e = starts_ends[i], starts_ends[i+1] 
        X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
        mus[i,:] = np.mean(X,1)    
    
    c =  np.corrcoef(mus)
    
    #zero out diagnals
    off_diag = []
    for j in xrange(c.shape[0]-1):
        off_diag.append(c[j,j+1])
    
    return c, off_diag

def get_correlated_segments(all_starts_ends, g, contig, r_cutoff):
    """
    take a set of coordinates representing starts and ends
    over a locus and cluster them into contiguous chunks that
    likely represent the same underlying call as a result of their
    high correlation
    """
    
    all_starts_ends = sorted(np.unique(all_starts_ends))
    
    c, off_diag = get_correlation_matrix(all_starts_ends, g, contig)
    #print all_starts_ends 
    original_c = c
    prev_gts = None
    while(np.amax(off_diag)>r_cutoff):
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
        c, off_diag = get_correlation_matrix(all_starts_ends, g, contig)

    s_e_tups = []
    for i in xrange(len(all_starts_ends)-1):
        s_e_tups.append([all_starts_ends[i], all_starts_ends[i+1]])
    
    return s_e_tups, original_c


def merge_correlated_calls(overlapping_call_clusts, g, contig, r_cutoff = 0.65):
    
    """
    take overlapping call clusts and determine the correlated segments
    if a whole block is correlated, it represents a single locus call
    """
        
    all_starts_ends = []
    min_s, max_e = 9e9, -1
    for clust in overlapping_call_clusts: 
        s,e = clust.get_med_start_end()
        min_s = min(min_s, s)
        max_e = max(max_e, e)
        all_starts_ends.append(s)
        all_starts_ends.append(e)
    
    s_e_segs, c = get_correlated_segments(all_starts_ends, g, contig, r_cutoff)
    print s_e_segs
    """
    1. remove crappy segs
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
        if gX.n_clusts > 1:
            if np.all(prev_labels == gX.labels):
                CNV_segs[-1][1] = e 
            else:
                CNV_segs.append([s,e])
                CNV_gXs.append(gX)
          
        prev_labels = gX.labels
    
    """
    remove calls (from all the clusts) where there are no CNV LOCI   
    v_calls = []
    for clust in overlapping_call_clusts:
        s, e =  clust.get_med_start_end()
        if overlap(s, e, CNV_segs):
            v_calls.append(clust)
    """
   
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

    """
    OK, now, get all the individual calls out of this and quickly recluster
    
    now - all these calls are the genotypes, so output
    region 1 - gts of indivs w/ that call, and then others, obviously don't have it
    region 2 - ditto 
    
    """
    
    """
    return a set of calls
    """

    if len(CNV_segs) <= 1 or non_adjacent(CNV_segs):
        return s_e_segs, CNV_segs, cnv_segs_by_indiv, c, v_calls, True
    else:
        for seg in CNV_segs:
            s, e = seg
            X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
            gX = g.GMM_genotype(X)
            Xs, s_idx_s, s_idx_e = g.get_sunk_gt_matrix(contig, s, e)
            gXs = g.GMM_genotype(Xs)
            g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, fn="./test/%d_%d_x.pdf"%(s+1, e))
        return s_e_segs, CNV_segs, cnv_segs_by_indiv, c, v_calls, False
    
def overlap(s, e, segs):
    for s_e in segs: 
        s_start, s_end = s_e
        if (s<=s_end and s>=s_start) or (e>=s_start and e<=s_end):
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
    g.plot(gX, gXs, contig, s, e, idx_s, idx_e, s_idx_s, s_idx_e, fn="./test/%d_%d.pdf"%(s, e))
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
        
        self.indivs = indivs 
        self.best_idx = self.bics.index(self.min_bic) 
        self.n_clusts = self.params[self.best_idx]
    
    def get_gts_by_indiv(self):
        
        cp_2_thresh=1.0
        m = mode(self.labels)[0]
        
        label_to_mu = {}
        mu_to_labels = {}

        all_labels = []
        all_mus = []
        
        for l in np.unique(self.labels):
            
            label_to_mu[l] = np.mean(self.mus[self.labels==l])
            mu_to_labels[np.mean(self.mus[self.labels==l])] = l
            all_labels.append(l)
            all_mus.append(np.mean(self.mus[self.labels==l]))
            
        all_labels = np.array(all_labels)
        all_mus = np.array(all_mus)
        
        mu_args = np.argsort(all_mus) 
          
        ordered_labels = all_labels[mu_args]
        ordered_mus = all_mus[mu_args] 
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
            new_labels_to_gt = {}
            for l, gt in labels_to_gt.iteritems():
                new_labels_to_gt[l] = gt+1
            labels_to_gt = new_labels_to_gt
        
        # finally double check the cps make sense...
        #sorted_gts = np.array(sorted(labels_to_gt.values()))
        #
        #l_eval = lambda x: np.sum(np.pow(x[0] - x[1]),2)
        #
        #v = l_eval([sorted_gts, ordered_mus])

        #for i xrange(1,len(sorted_gts)):
        #    if l_eval([sorted_gts[i:]+1,ordered_mus]) < v:
        #        labels_to_gt
        #        sorted_gts = np.array(sorted(labels_to_gt.values()))
        #        v = l_eval([sorted_gts, ordered_mus])
        gts_by_indiv = {}
        for i, indiv in enumerate(self.indivs):  
            gts_by_indiv[indiv] = labels_to_gt[self.labels[i]] 
        
        return gts_by_indiv
            
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
        
def output(call_clust, g, contig, s, e, F_gt, F_call, include_indivs=None, plot=False):

    X, idx_s, idx_e = g.get_gt_matrix(contig, s, e)
    gX = g.GMM_genotype(X, include_indivs)
    if gX.n_clusts == 1:
        return

    F_call.write("%s\t%d\t%d\n"%(contig, s, e))
    g.output(F_gt, gX, contig, s, e)


class genotyper:
    
    def setup_output(self, FOUT):
        outstr = "contig\tstart\tend\t%s\n"%("\t".join(self.indivs))
        FOUT.write(outstr)

    def output(self, FOUT, gX, contig, s, e):
        
        outstr = "%s\t%d\t%d"%(contig, s, e)

        gts_by_indiv = gX.get_gts_by_indiv()
        
        ordered_gts = []
        for indiv in self.indivs:
            if indiv in gts_by_indiv:
                ordered_gts.append(gts_by_indiv[indiv])
            else:
                ordered_gts.append(-1)
        
        outstr = "%s\t%s\n"%(outstr, "\t".join("%d"%gt for gt in ordered_gts))
        FOUT.write(outstr)
    
    def init_on_indiv_DTS_files(self, **kwargs):

        g = gglob.init_from_DTS(**kwags)

        self.indivs = g.indivs
        self.wnd_starts = g.wnd_starts
        self.wnd_ends = g.wnd_ends
        self.sunk_wnd_starts = g.sunk_wnd_starts
        self.sunk_wnd_ends = g.sunk_wnd_ends
        
        self.cp_matrix = g.cp_matrix
        self.sunk_cp_matrix = g.sunk_cp_matrix

    
    def init_on_gglob(self, contig, fn_gglob):
        
        g = gglob.init_from_gglob_dir(contig, fn_gglob)
        
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
        
        self.contig = contig
        self.indivs = []
        self.wnd_starts = None
        self.wnd_ends = None
        self.sunk_wnd_starts = None
        self.sunk_wnd_ends = None
        
        self.cp_matrix = None
        self.sunk_cp_matrix = None
        
        if self.gglob_dir:
            self.init_on_gglob(self.gglob_dir, self.contig) 
        else:
            self.init_on_indiv_DTS_files(self, **kwargs)

        print >>stderr, "loading genomes..."
        t = time.time()
        
        print >>stderr, "done (%fs)"%(time.time()-t)
       
    def addGMM(self, gmm, ax, X):
        
        G_x=np.arange(0,max(X)+1,.1)
        l = gmm.means.shape[0] 
        print l
        for i in xrange(l):
            c = cm.hsv(float(i)/l,1)
            mu = gmm.means[i,0]
            var = gmm.covars[i][0][0]
            print mu, var

            G_y = mlab.normpdf(G_x, mu, var**.5)*gmm.weights[i]
            ax.plot(G_x,G_y,color=c)
            ax.plot(mu,-.001,"^",ms=10,alpha=.7,color=c)
            
    def aug_dendrogram(self, ax, ddata):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            ax.plot(x, y, 'ro')
            ax.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')

    def plot(self, gX, gXs, chr, start, end, idx_s, idx_e, s_idx_s, s_idx_e, fn="test_gt.pdf"):
        
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
        
        n, bins, patches = axarr[1,1].hist(cps,alpha=.9,ec='none',normed=1,color='r',bins=len(cps)/20)
        self.addGMM(gX.gmm, axarr[1,1], cps)
        fig.sca(axarr[0,2]) 
        dendro = hclust.dendrogram(gX.Z, orientation='right')
        #self.aug_dendrogram(axarr[0,2], dendro)
        
        n, bins, patches = axarr[1,0].hist(sunk_cps,alpha=.9,ec='none',normed=1,color='g',bins=len(cps)/20)
        self.addGMM(gXs.gmm, axarr[1,0], sunk_cps)
        axarr[0,1].plot(gX.params, gX.bics, 'ro-')
        axarr[0,1].plot(gXs.params, gXs.bics, 'go-')
        
        fig.sca(axarr[1,2]) 
        dendro = hclust.dendrogram(gXs.Z, orientation='right')
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
        
        print "BEGIN:", start_idx, end_idx 
        if end_idx<=start_idx:
            end_idx = start_idx+1
        elif end_idx-start_idx>1:
            start_idx+=1
        print "END:", start_idx, end_idx 

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
    
    def fit_GMM(self, X, init_means, init_vars, init_weights):
    
        n_components = len(init_means)
        gmm = mixture.GMM(n_components, 'spherical')
        gmm.means = np.reshape(np.array(init_means),(len(init_means),1))
        gmm.weights = np.array(init_weights)
        
        #vars = np.array([v[0][0] for v in gmm.covars])
        #gmm.covars = np.reshape()

        gmm.fit(X, n_iter=1000, init_params='c')
        labels = gmm.predict(X)
        
        bic = -2*gmm.score(X).sum() + (3*n_components)*np.log(X.shape[0])
        aic = -2*gmm.score(X).sum() + 2*(3*n_components)
        
        return gmm, labels, bic 

    def GMM_genotype(self, X, include_indivs = None, FOUT = None):
        """
        GMM genotyping
        begin by subsetting X to those indivs you want
        #cv_types = ['spherical', 'tied', 'diag', 'full']
        """
         
        if include_indivs:
            l = len_include_indivs
            new_X = np.zeros(l, X.shape[1])
            
            j = 0
            for i in xrange(X.shape[0]):
                if self.indivs[i] in include_indivs:
                    new_X[j] = X[i] 
                    j+=1
            X=new_X
                    
        mus = np.mean(X,1)
        mus = np.reshape(mus, (mus.shape[0],1))
        
        dist_mat = dist.pdist(mus)
        Z = hclust.linkage(mus, method='centroid', metric='euclidean')
        params, bics, gmms, all_labels = [], [], [], []
        
        prev_grps = np.array([])
        for k in np.arange(.2, 0.7,  0.01):
            grps = hclust.fcluster(Z, k, criterion='distance')
            if np.all(grps == prev_grps): continue

            init_mus, init_vars, init_weights = self.initialize(mus, grps) 

            gmm, labels, ic = self.fit_GMM(mus, init_mus, init_vars, init_weights)

            params.append(len(init_mus))
            bics.append(ic)
            gmms.append(gmm)
            all_labels.append(labels)
            prev_grps = grps 

        grps = np.zeros(mus.shape[0])
        init_mus, init_vars, init_weights = self.initialize(mus, grps) 
        gmm, labels, ic = self.fit_GMM(mus, init_mus, init_vars, init_weights)
        params.append(len(init_mus))
        bics.append(ic)
        gmms.append(gmm)
        all_labels.append(labels)
        
        #overlaps = self.pw_GMM_overlap(gmm)
        
        #if FOUT: 
        #    FOUT.write(" ".join(["%f"%(max(o)) for o in overlaps]))
        idx = np.argmin(bics)
        gmm = gmms[idx]
        labels = all_labels[idx]
        if include_indivs == None: 
            include_indivs = self.indivs
            
        return GMM_gt(X, gmm, labels, Z, params, bics, include_indivs)
        #return gmm, labels, Z, [params, bics]
        #return gmm, labels, Z, [[len(init_mus)], [ic]]
    
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
        
    def pw_GMM_overlap(self, gmm):
        
        overlaps = []       
        mus = gmm.means[:,0]
        vars = np.array([v[0][0] for v in gmm.covars])
        weights = np.array(gmm.weights)
        order = np.argsort(mus)
        mus = mus[order]
        weights = weights[order]
        vars = vars[order]
        
        for i in xrange(mus.shape[0]-1):
            mu1 = mus[i]
            mu2 = mus[i+1]
            v1 = vars[i]
            v2 = vars[i+1]
            sd_max = np.sqrt(max(v1, v2))
            mn = min(mu1, mu2) - 10*sd_max
            mx = max(mu1, mu2) + 10*sd_max
            xs = np.arange(mn,mx,0.01)
            o = np.sum(np.min(np.c_[norm.pdf(xs,loc=mu1,scale=v1)*weights[i], 
                                    norm.pdf(xs,loc=mu2,scale=v2)*weights[i+1]], 1))  * 0.01
            overlaps.append([o,o/weights[i],o/weights[i+1]])
        return overlaps
    
    def __plot(self, Z, cutoff):
        
        dendro = hclust.dendrogram(Z, orientation='right', labels = lbls, color_threshold = cutoff)
        grps = hclust.fcluster(Z, cutoff, criterion='distance')
        
        plt.gcf().set_size_inches(14,6)
        ax=plt.gcf().gca()
        ax.set_position([.05,0.05,.3,.9])

        ax2 = plt.gcf().add_axes([.55,.05,.4,.9])
        k = 0
        colors = ['b','g','r','c','m','y','k']
        n_colors = len(colors)

        for clust in Set(list(grps)):
            for idx in np.where(grps == clust)[0]:
                ax2.plot([self.calls[idx]['start'], self.calls[idx]['end']],
                         [k,k], 
                         lw=1,
                         color=colors[(clust-1)%n_colors])
                k+=1.2
          
        ax2.set_xlim([min(self.all_starts),max(self.all_ends)])
        ax2.set_ylim([-1,k+1])

        plt.savefig('test.png')
        plt.gcf().clear()

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
    
