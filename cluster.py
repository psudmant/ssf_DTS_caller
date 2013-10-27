"""
classes and helper functions for clustering calls
"""
import pandas as pd
from sys import stderr
import copy
from collections import defaultdict
import numpy as np
import random
import scipy.cluster.hierarchy as hclust
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist 
from sets import Set
import networkx as nx
import time

class callset_table:
    """
    a callset table must have columns with a test, reference
        indiv_ref
        indiv_test
        *chr
        *start
        *end
        mu
        p
        window_size
    """

    def __init__(self, fn_table):
        print >>stderr, "loading table..."
        t=time.time()
        self.pd_table =  pd.read_csv(fn_table, 
                                     sep = '\t', 
                                     header=0, 
                                     compression = 'gzip')
        print >>stderr, "done (%fs)"%(time.time()-t)

    def filter(self, p_cutoff, size_cutoff):
        """
        apply size cutoff (no calls of size 1 window!)
        """
        print >>stderr, "parsing calls with p<=%f and l>=%d..."%(p_cutoff, size_cutoff)
        self.pd_table = self.pd_table[(self.pd_table['p']<p_cutoff)]
        self.pd_table = self.pd_table[(self.pd_table['window_size']>=size_cutoff)]
        #self.pd_table = self.pd_table[( (self.pd_table['window_size']>=size_cutoff) |
        #                                (self.pd_table['p']>0) )]
        print >>stderr, "done"
    
    def sort(self):
        """
        sort on chr, start, end
        """
        print >>stderr, "sorting..."
        t=time.time()
        self.pd_table = self.pd_table.sort(['chr', 'start', 'end', 'p'])
        print >>stderr, "done (%fs)"%(time.time()-t)

class cluster_calls:
    
    def __init__(self, callset_table):
        """
        callset_table is a pandas
        """
        self.callset_table = copy.deepcopy(callset_table)
        self.callset_table.sort()
        self.overlapping_calls_by_chr, self.n_overlapping = self.get_overlapping_calls()
        print >>stderr, "%d overlapping calls"%self.n_overlapping
    
    def resolve_overlapping_clusters(self, verbose=False):
        """
        resolve overlapping clusters into individual calls 
        let the calls have likelihoods
        """
        ll_cutoff = -6.0
        print >>stderr, "resolving breakpoints..."
        final_calls = []
        
        for chr, overlapping_calls in self.overlapping_calls_by_chr.iteritems():
            for overlap_cluster in overlapping_calls:
                resolved_calls = overlap_cluster.resolve(0.85, ll_cutoff, min_size=2)
                final_calls += resolved_calls
                if verbose:
                    for call in resolved_calls:
                        call.print_verbose()
                    
        print >>stderr, "done"
        print >>stderr, "%d calls with likelihood <%f"%(len(final_calls), ll_cutoff)
        return final_calls


    def output_overlap_clusters(self, fn, id):
        """
        track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2
        itemRgb="On"
        chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
        """
        with open(fn,'w') as F:    
            print >>F, """track name="%s" description="%s" visibility=2 itemRgb="On" """%(id, id)
            for chr, overlapping_calls in self.overlapping_calls_by_chr.iteritems():
                for overlap_cluster in overlapping_calls:
                    rgb_str = "%d,%d,%d"%(tuple([random.randrange(0,255,1) for r in xrange(3)]))
                    for call in overlap_cluster.calls:
                        print >>F, "%s\t%d\t%d\t%f\t%f\t+\t%d\t%d\t%s"%(overlap_cluster.chr, 
                                call['start'],
                                call['end'],
                                call['p'],
                                call['p'],
                                call['start'],
                                call['end'],
                                rgb_str)
    
    def get_overlapping_calls(self):
        """
        glob up all those calls that physically overlap
        assumes self.callset_table is sorted
        """
        curr_max=self.callset_table.pd_table.iloc[0]['end']
        curr_chr=self.callset_table.pd_table.iloc[0]['chr']
        
        overlapping_calls_by_chr = defaultdict(list)
        #overlapping_calls = []
        
        curr_cluster = call_cluster(curr_chr) 

        overlapping_calls_by_chr[curr_chr].append(curr_cluster)
        n=1

        for idx, row in self.callset_table.pd_table.iterrows():
            chr, start, end = row['chr'], row['start'], row['end']
            if start<curr_max and curr_chr == chr:
                curr_cluster.add(row)
            else:
                curr_cluster = call_cluster(chr) 
                curr_chr = chr
                curr_max = end
                curr_cluster.add(row)
                overlapping_calls_by_chr[curr_chr].append(curr_cluster)
                n+=1
            
            curr_max = max(end, curr_max) 
            curr_chr = chr 
        
        return overlapping_calls_by_chr, n


def intersect(si,ei,sj,ej):
    
    if  si < ej and sj < ei: 
        return True
    else:
        return False

def get_overlapping_call_clusts(recip_overlap_clusts):
    
    G = nx.Graph()
    
    for clust in recip_overlap_clusts:
        G.add_node(clust)

    for clust1 in recip_overlap_clusts:
        si, ei = clust1.get_best_start(), clust1.get_best_end()
        for clust2 in recip_overlap_clusts:
            sj, ej = clust2.get_best_start(), clust2.get_best_end()
            if intersect(si,ei,sj,ej):
                G.add_edge(clust1, clust2)
        
    return nx.connected_components(G)

class call_cluster:
    def __init__(self, chr):
        self.calls = []
        self.all_starts = []
        self.all_ends = []
        self.chr = chr
        self.size = 0
        self.log_likelihood =  None
        """
        'complex' refers to calls that have been 
        clustered (overlap) together, but reach significance
        self.compound = False 
        """
         
    def add(self, call):
        self.calls.append(call)
        self.all_starts.append(call['start'])
        self.all_ends.append(call['end'])
        self.size+=1
    def get_indiv_calls_str(self):
        ret_str = ""
        for call in self.calls:
            ret_str+="%s\n"%("\t".join([str(s) for s in [call["indiv_ref"], 
                                               call["indiv_test"], 
                                               call["chr"],
                                               call["start"], 
                                               call["end"], 
                                               call["mu"],
                                               call["p"], 
                                               call["window_size"]]]))
        return ret_str
            
    #def get_start(self):
    #    if size==1:
    #        return self.calls[0]['start']
    #    elif not selfcomplex

    def get_best_start(self):
        dd_starts = defaultdict(int)
        for call in self.calls:
            dd_starts[call['start']]+=np.log10(call['p'])
        
        return sorted(dd_starts.iteritems(), key=lambda x: x[1])[0][0]
    
    def get_best_end(self):
        dd_starts = defaultdict(int)
        for call in self.calls:
            dd_starts[call['end']]+=np.log10(call['p'])
        
        return sorted(dd_starts.iteritems(), key=lambda x: x[1])[0][0]

    #def get_med_start(self):
    #    return np.median(np.array(self.all_starts))
    #
    #def get_med_end(self):
    #    return np.median(np.array(self.all_ends))

    def get_log_likelihood(self, force=False):

        if self.log_likelihood == None or force:
            self.log_likelihood= np.sum(np.log10(np.array([call['p'] for call in self.calls])))
        return self.log_likelihood

    def get_range_str(self):
        return "%s:%d-%d"%(self.chr,min(self.all_starts), max(self.all_ends))
    
    def print_out(self):
        for c in self.calls:
            print "\t", c['chr'], c['start'], c['end'], np.log10(c['p']), c['window_size']

    def frac_overlap(self, i, j):
        si, ei = self.calls[i]['start'],  self.calls[i]['end']  
        sj, ej = self.calls[j]['start'],  self.calls[j]['end']  

        if not (si < ej and sj < ei): 
            return 0   

        overlap = min(ei,ej) - max(si, sj)
        return min(float(overlap)/(ei-si), 
                   float(overlap)/(ej-sj))
    
    def plot(self, Z, cutoff):
        print self.chr, min(self.all_starts), max(self.all_ends)
        lbls = ["%s:%d-%d"%(c['chr'], c['start'], c['end']) for c in self.calls]
        
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
        
    def cluster_by_recip_overlap(self, cutoff, ll_cutoff=0, plot=False):

        mat = np.zeros((self.size,self.size))
        for i in xrange(self.size):
            for j in xrange(self.size):
                #mat[i,j] = 1-((self.frac_overlap(i,j)+self.frac_overlap(j,i))/2.0)
                mat[i,j] = 1-(self.frac_overlap(i,j))

        
        mat=1.-dist.squareform( (1-mat) * (np.diag(np.repeat(1,self.size),0) -1) * -1)
        Z = hclust.linkage(mat, method='complete')
        grps = hclust.fcluster(Z, cutoff, criterion='distance')
        
        if plot:
            self.plot(Z, cutoff)
            raw_input()
        
        new_groups = []
        for grp in np.unique(grps):
            new_group = call_cluster(self.chr) 
            for idx in np.where(grps==grp)[0]:
                new_group.add(self.calls[idx])
            if new_group.get_log_likelihood()<ll_cutoff: new_groups.append(new_group)
        return new_groups
        
    def resolve(self, overlap_cutoff, ll_cutoff, min_size=1, do_plot=False):
        """
        though a bunch of calls may overlap, they may not be all 'good'
        or all the same call. Thus, further cluster overlapping calls
        into those with reciprocal overlap properties
        """
         
        if self.size <min_size: 
            return []

        #single event (in the event that we ARE allowing them, by default, not (see above catch)     
        if self.size == 1:
            return self.get_log_likelihood < ll_cutoff and [CNV_call([self],self.chr)] or []
                                                           
        recip_overlap_clusts = self.cluster_by_recip_overlap(overlap_cutoff, ll_cutoff=ll_cutoff, plot=do_plot)
        
        
        overlapping_call_clusts = get_overlapping_call_clusts(recip_overlap_clusts)
        CNV_calls = [CNV_call(c, self.chr) for c in overlapping_call_clusts]
        
        if do_plot:
            self.print_out()
            for call_clust in sorted(recip_overlap_clusts, key=lambda x: x.get_log_likelihood()):
                print "\t", call_clust.get_range_str(), call_clust.size, call_clust.get_log_likelihood()
            for call in CNV_calls:
                print ">final call:", call.print_str()
        
        return  CNV_calls

class CNV_call:
    """
    the purpose of this is, after filtering clustered calls
    we want to have a final object represneting a CNV that's been called
    It may (usually) will consist of ONE call_cluster, however, 
    sometimes, multiple call_clusters overlap (just not enough to be clustered 
    by recip overlap). In this case, the CNV call holds all of these, 
    and the coordinates are the min and max of the med_start/ends
    """
    def __init__(self, clustered_calls, chr):
        self.clustered_calls = clustered_calls
        self.chr = chr
        
        best_ll = 1
        for clust in clustered_calls:
            if clust.get_log_likelihood() < best_ll:
                self.best_start = clust.get_best_start()
                self.best_end = clust.get_best_end()
        
        self.best_log_likelihood=best_ll
    
    def print_str(self):
        
        outstr=""
        for clust in self.clustered_calls:
            start = clust.get_best_start()
            end = clust.get_best_end()
            outstr+="%s\t%d\t%d\t%f\n"%(self.chr, start, end, clust.get_log_likelihood())
        return outstr
    
    def get_indiv_calls_str(self):
        
        outstr=""
        for clust in self.clustered_calls:
            outstr+=clust.get_indiv_calls_str()
        return outstr
    
    def print_verbose(self):
        print "%s\t%d\t%d\tll:%f\t%d clusts"%(self.chr, 
                                              self.start, 
                                              self.end, 
                                              self.log_likelihood, 
                                              len(self.clustered_calls))
        for clust in self.clustered_calls:
            print clust.chr, clust.get_best_start(), clust.get_best_end()
            clust.print_out()


    """ 
    def assess(self):
        self.calc_combined_p()
        self.get_simple_bounds()
       
    def get_simple_bounds(self):
        self.mode_start = int(stats.mode(np.array(self.all_starts))[0][0])
        self.mode_end = int(stats.mode(np.array(self.all_ends))[0][0])
        if self.mode_start > self.mode_end:
            print "||||||||WARNING START > END||||||||||||"
            print self.chr, self.mode_end, self.mode_start, self.mode_start - self.mode_end
            self.print_out()
            e = self.mode_end
            self.mode_end = self.mode_start
            self.mode_start = e

    def calc_combined_p(self):
        self.combined_log_p = 0 
        for c in self.calls:
            if c['p']==0: c['p']=1e-20
            self.combined_log_p += -math.log(c['p'])
    """ 





"""
def cluster_calls(calls):
    
    curr_max=calls.iloc[0]['end']
    curr_chr=calls.iloc[0]['chr']
    
    call_clusters_by_chr = defaultdict(list)
    curr_cluster = call_cluster() 

    call_clusters_by_chr[curr_chr].append(curr_cluster)
    n=1

    for idx, row in calls.iterrows():
        chr, start, end = row['chr'], row['start'], row['end']
        
        if start<curr_max and curr_chr == chr: #and curr_cluster.overlap(start,end)>0.4:
            curr_cluster.add(row)
        else:
            curr_cluster.assess()
            curr_cluster = call_cluster() 
            curr_cluster.add(row)
            call_clusters_by_chr[curr_chr].append(curr_cluster)
            n+=1
        
        curr_max = end 
        curr_chr = chr 
    curr_cluster.assess()
    
    return call_clusters_by_chr, n


def output_bed_of_clusters(call_clusters_by_chr, fn_out):
    
    F = open(fn_out,'w')
        
    for chr, clusters in call_clusters_by_chr.iteritems():
        print chr
        for cluster in sorted(clusters,key=lambda x: -x.combined_log_p):
            print >>F, "%s\t%d\t%d\t%f"%(cluster.chr, cluster.mode_start, cluster.mode_end, cluster.combined_log_p) 
"""






