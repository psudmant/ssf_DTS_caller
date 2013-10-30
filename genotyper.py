import glob
from sys import stderr
import numpy as np

from wnd_cp_data import wnd_cp_indiv

import time
import matplotlib.pyplot as plt
import matplotlib.colors as mCols
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

class genotyper:

    def __init__(self, 
                 dts_dir, 
                 sunk_dts_dir, 
                 fn_contigs, 
                 fn_sunk_contigs, 
                 wnd_size,
                 contig,
                 plot_dir,
                 F_fnToIndiv=lambda x: x.split("/")[-1].replace("500_bp_",""),
                 limit_to_n=None):
        
        self.i_by_indiv = {}
        self.indivs = []
        self.wnd_starts = None
        self.wnd_ends = None
        self.sunk_wnd_starts = None
        self.sunk_wnd_ends = None
        
        self.cp_matrix = None
        self.sunk_cp_matrix = None
        
        self.has_sunk_cps = False
        
        if fn_sunk_contigs != None:
            self.has_sunk_cps = True

        self.contig = contig
        self.plot_dir = plot_dir
        
        DTS_prefix="500_bp_" 

        print >>stderr, "loading genomes..."
        t = time.time()
        
        fn_DTSs = glob.glob("%s/%s*"%(dts_dir, DTS_prefix))
        
        n_indivs = limit_to_n and limit_to_n or len(fn_DTSs)
        
        rand_wnd_cp = wnd_cp_indiv(fn_DTSs[0], fn_contigs, wnd_size)
        self.wnd_starts, self.wnd_ends = rand_wnd_cp.get_wnds_by_chr(contig)
        self.cp_matrix = np.zeros((n_indivs, self.wnd_starts.shape[0]))

        if self.has_sunk_cps:
            fn_sunk_DTSs = glob.glob("%s/%s*"%(sunk_dts_dir, DTS_prefix))
            rand_sunk_wnd_cp = wnd_cp_indiv(fn_sunk_DTSs[0], fn_sunk_contigs, wnd_size)
            self.sunk_wnd_starts, self.sunk_wnd_ends = rand_sunk_wnd_cp.get_wnds_by_chr(contig)
            self.sunk_cp_matrix = np.zeros((n_indivs, self.sunk_wnd_starts.shape[0]))

        if 0:
            for i, fn_DTS in enumerate(fn_DTSs):
                if limit_to_n and i>=limit_to_n: break

                indiv = F_fnToIndiv(fn_DTS)
                wnd_cp = wnd_cp_indiv(fn_DTS,
                                      fn_contigs,
                                      wnd_size)
                
                self.indivs.append(indiv)
                self.cp_matrix[i,:] = wnd_cp.get_cps_by_chr(contig,correct=True) 

                if self.has_sunk_cps:
                    wnd_cp = wnd_cp_indiv("%s/%s%s"%(sunk_dts_dir, DTS_prefix, indiv),
                                          fn_sunk_contigs,
                                          wnd_size)
                    self.sunk_cp_matrix[i,:] = wnd_cp.get_cps_by_chr(contig, correct=True) 
                                                                         
            np.save("./tmp_cp_matrix", self.cp_matrix)
            np.save("./tmp_sunk_cp_matrix", self.sunk_cp_matrix)
        else:
            self.cp_matrix = np.load("./tmp_cp_matrix.npy") 
            self.sunk_cp_matrix = np.load("./tmp_sunk_cp_matrix.npy") 

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
            
    def plot(self, X, Xs, gmmX, gmmXs, v_bndsX, v_bndsXs, chr, start, end):
        
        cps = np.mean(X, 1)
        sunk_cps = np.mean(Xs, 1)
        
        plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
        fig, axarr = plt.subplots(2, 2)
        fig.set_figwidth(11)
        fig.set_figheight(8.5)
        axescolor  = '#f6f6f6'
        #h_margin = 0.05
        #v_margin = 0.05
        #plot_height=1-2*v_margin
        #plot_width=((1-2*h_margin)-2*h_margin)/3.0
        #bottom=1-h_margin-plot_height
        #plot_rects = []
        #plot_rect=[h_margin,bottom,plot_width,plot_height]
        #ax = fig.add_axes(plot_rect)
        #
        #ax.plot(cps, sunk_cps, 'ro', alpha=0.5)
        
        axarr[0,0].plot(cps, sunk_cps, 'ro', alpha=0.2)
        axarr[0,0].set_xlim(-0.10,max(cps)+1)
        axarr[0,0].set_ylim(-0.10,max(sunk_cps)+1)
        
        n, bins, patches = axarr[0,1].hist(cps,alpha=.9,ec='none',normed=1,color='r',bins=len(cps)/20)
        self.addGMM(gmmX, axarr[0,1], cps)
        n, bins, patches = axarr[1,0].hist(sunk_cps,alpha=.9,ec='none',normed=1,color='g',bins=len(cps)/20)
        self.addGMM(gmmXs, axarr[1,0], sunk_cps)
        axarr[1,1].plot(v_bndsX[0], v_bndsX[1], 'ro-')
        axarr[1,1].plot(v_bndsXs[0], v_bndsXs[1], 'go-')
        
        fig.savefig("%s/%s-%d-%d.png"%(self.plot_dir, chr, start, end))
        plt.close()


    def _plot(self, X, labels, chr, start, end):

        cps = np.mean(X, 1)
        print cps.shape
        plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
        fig = plt.figure()
        fig.set_figwidth(8)
        fig.set_figheight(6)
        axescolor  = '#f6f6f6'
        h_margin = 0.05
        v_margin = 0.05
        plot_height=1-2*v_margin
        plot_width=1-2*h_margin
        bottom=1-h_margin-plot_height
        plot_rect=[h_margin,bottom,plot_width,plot_height]

        hist_ax = fig.add_axes(plot_rect)
        n, bins, patches = hist_ax.hist(cps,alpha=.9,ec='none',normed=1,bins=len(cps)/20)
        
        uniq_labels = list(Set(labels))
        n_labels = len(uniq_labels)
        
        G_x=np.arange(0,max(cps)+2,.1)
        
        for i, label in enumerate(uniq_labels):

            c=cm.hsv(float(i)/n_labels,1)
            mu = np.mean(cps[labels==label])
            var = np.var(cps[labels==label])
            frac = np.sum(labels==label)/float(len(labels)) 
            n_in_label = np.sum(labels==label)
            G_y = mlab.normpdf(G_x, mu, var**.5)*frac
            hist_ax.plot(G_x,G_y,color=c)
            hist_ax.plot(mu,-.001,"^",ms=10,alpha=.7,color=c)
            print "%d cluster at %f"%(n_in_label, mu) 
        
        fig.savefig("%s/%s-%d-%d.png"%(self.plot_dir, chr, start, end))
        plt.close()
        
    def get_gt_matrix(self, contig, start, end, vb=False):
        assert contig == self.contig
        start_idx = np.searchsorted(self.wnd_starts, start)
        end_idx = np.searchsorted(self.wnd_ends, end)
        
        X = self.cp_matrix[:,start_idx:end_idx]
        if vb:
            print "idxs:", start_idx, end_idx
            print X
        return X
    
    def get_sunk_gt_matrix(self, contig, start, end, vb=False):
        assert contig == self.contig
        assert self.has_sunk_cps == True

        start_idx = np.searchsorted(self.sunk_wnd_starts, start)
        end_idx = np.searchsorted(self.sunk_wnd_ends, end)
        
        X = self.sunk_cp_matrix[:,start_idx:end_idx]
        if vb:
            print "idxs:", start_idx, end_idx
            print X
        return X
    
    def get_gt_matrix_mu(self, contig, start, end):
        assert contig == self.contig
        start_idx = np.searchsorted(self.wnd_starts, start)
        end_idx = np.searchsorted(self.wnd_ends, end)
        X = self.cp_matrix[:,start_idx:end_idx]
        
        return X

    def s_score(self, X, labels):
        return metrics.silhouette_score(X, labels) 
    

    def fit_GMM(self, X, init_means, init_weights):
    
        n_components = len(init_means)
        gmm = mixture.GMM(n_components, 'spherical')
        gmm.means = np.reshape(np.array(init_means),(len(init_means),1))
        gmm.weights = np.array(init_weights)
        gmm.fit(X, n_iter=1000, init_params='c')
        labels = gmm.predict(X)
        
        bic = -2*gmm.score(X).sum() + (3*n_components)*np.log(X.shape[0])
        aic = -2*gmm.score(X).sum() + 2*(3*n_components)
        
        return gmm, labels, bic 

    def GMM_genotype(self, X):
        """
        GMM genotyping
        #cv_types = ['spherical', 'tied', 'diag', 'full']
        """
        mus = np.mean(X,1)
        #mus = np.median(X,1)
        mus = np.reshape(mus, (mus.shape[0],1))
        
        round_mus = np.around(mus) 
        init_mus = list(np.unique(round_mus))
        
        init_weights = []
        for comp in init_mus:
            init_weights.append(float(np.sum(round_mus==comp))/len(round_mus))
        
        gmm, labels, ic = self.fit_GMM(mus, init_mus, init_weights)
        
        print self.pw_GMM_overlap(gmm)
        return gmm, labels, [[len(init_mus)], [ic]]
    
    def pw_GMM_overlap(self, gmm):
        
        overlaps = []       
        mus = gmm.means[:,0]
        vars = np.array([v[0][0] for v in gmm.covars])
        order = np.argsort(mus)
        mus = mus[order]
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
            o = np.sum(np.min(np.c_[norm.pdf(xs,loc=mu1,scale=v1), 
                                    norm.pdf(xs,loc=mu2,scale=v2)], 1))  * 0.01
            print mn, mx, mu1, mu2, v1, v2, o
            overlaps.append(o)
        return overlaps

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
    
