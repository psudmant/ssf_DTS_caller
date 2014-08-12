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
import matplotlib
import scipy.ndimage as ndi
import scipy.spatial.distance as dist 
from sets import Set
import networkx as nx
import time
import pysam
import math
import pdb

from sets import Set

from genotyper import genotyper

class callset_table(object):
    """
        *chr
        *start
        *end
        [indiv_ref
        indiv_test
        mu
        p
        window_size]
    """

    def __init__(self, fn_table):
        print >>stderr, "loading table..."
        t=time.time()
        self.pd_table =  pd.read_csv(fn_table, 
                                     sep = '\t', 
                                     header=0, 
                                     compression = 'gzip')
        print >>stderr, "done (%fs)"%(time.time()-t)
    
    def filter_by_chr(self, contig):

        print >>stderr, "filtering calls to only %s"%(contig)
        self.pd_table = self.pd_table[self.pd_table['chr']==contig]

    def filter_by_gsize(self, max_size):

        print >>stderr, "filtering calls >%dbp"%(max_size)
        self.pd_table = self.pd_table[(self.pd_table['end']-self.pd_table['start'])<=max_size]
    
    def sort(self):
        """
        sort on chr, start, end
        """
        print >>stderr, "sorting..."
        t=time.time()
        self.pd_table = self.pd_table.sort(['chr', 'start', 'end'])
        print >>stderr, "done (%fs)"%(time.time()-t)

class simple_callset_table(callset_table):

    def __init__(self, fn_table):
        super(simple_callset_table, self).__init__(fn_table)



class indiv_callset_table(callset_table):

    def __init__(self, fn_table):
        super(indiv_callset_table, self).__init__(fn_table)
        
    def output(self, fn_out):
        print >>stderr, "outputting bed file of filtered calls..."
        indiv=fn_out.split("/")[-1].split(".")[0]
        Fout = open(fn_out,'w')
        Fout.write("""track name="%s_indiv_calls" description="%s_indivs_calls" visibility=2 itemRgb="On"\n"""%(indiv, indiv))
        for idx, row in self.pd_table.iterrows():
            chr = row['chr']
            start = row['start']
            end = row['end']
            p = np.log(row['p'])
            ref = row['indiv_ref']
            mu = row['mu']
            Fout.write("%s\t%d\t%d\t%f_%s\t%f\t+\t%d\t%d\t%s\n"%(chr,start,end,p,ref,mu,start,end,"0,0,0"))
        Fout.close()

    def filter(self, p_cutoff, size_cutoff, single_window_cutoff, divide_by_mu=False, min_mu=0.5):
        """
        apply size cutoff (no calls of size 1 window!)
        """
        print >>stderr, "parsing calls with p<=%.4g and l>=%d..."%(p_cutoff, size_cutoff)
        print >>stderr, "initial table size: %d"%(self.pd_table.shape[0])
        if divide_by_mu: 
            where_p_sig=((self.pd_table['p']/np.exp((self.pd_table['mu']**2)))<p_cutoff)
            where_mu_sig=(np.absolute(self.pd_table['mu'])>min_mu)
            self.pd_table = self.pd_table[(where_p_sig | where_mu_sig)]
            print "p-sig:", np.sum(where_p_sig)
            print "mu-sig:", np.sum(where_mu_sig)
        else:
            self.pd_table = self.pd_table[(self.pd_table['p']<p_cutoff)]

        self.pd_table = self.pd_table[(self.pd_table['window_size']>=size_cutoff)]
        
        ##for windows of size 1, if they exist, ensure they exceed the max_mu
        #self.pd_table = self.pd_table[((self.pd_table['window_size']>1)|(np.absolute(self.pd_table['mu'])>=single_window_cutoff))]
        #self.pd_table = self.pd_table[( (self.pd_table['window_size']>=size_cutoff) |
        #                                (self.pd_table['p']>0) )]
        print >>stderr, "final table size: %d"%(self.pd_table.shape[0])
        print >>stderr, "done"

class simple_call:
        
    def __init__(self, contig, start, end, ref_indivs, resolved_calls):
        self.contig = contig
        self.start = start
        self.end = end
        self.ref_indivs  = ref_indivs #indivs CALLED over region
        self.resolved_calls = resolved_calls
        self.length = self.end-self.start
        self.total_ll = np.sum(np.array([call.get_log_likelihood() for call in self.resolved_calls]))


class indiv_cluster_calls:
    """
    a class to cluster calls made on one individual
    compared to multiple dCGH references
    """
    
    def __init__(self, callset_table):
        """
        callset_table is a pandas
        """
        self.callset_table = copy.deepcopy(callset_table)
        self.callset_table.sort()
        print >> stderr, "getting overlapping elements..."
        t=time.time()
        self.overlapping_calls_by_chr, self.n_overlapping = self.get_overlapping_calls()
        print >>stderr, "done (%fs)"%(time.time()-t)
        print >>stderr, "%d overlapping calls"%self.n_overlapping
    
    def get_curr_chr_dCGHs(self,chr, dCGHs):
        
        curr_dCGHs = {}
        for key,d in dCGHs.iteritems():
            curr_dCGHs[key] = d.get_cp_ratio_by_chr(chr)
        
        return curr_dCGHs


    def resolve_overlapping_clusters(self, ll_cutoff, 
                                           tbx_dups, 
                                           indiv_id, 
                                           indiv_DTS, 
                                           ref_DTSs, 
                                           dCGHs, 
                                           gglob_dir, 
                                           out_viz_dir,
                                           verbose=False, 
                                           min_overlapping=2, 
                                           subset_indivs=None,
                                           min_d=0):
        """
        resolve overlapping clusters into individual calls 
        let the calls have likelihoods
            1. do the recip overlap cluster - clusters very similar calls w/ similar break-points
            2. make sure there are at least 3 calls in there (w/ similar breakpoints)
            3. make sure those calls sum to a log likelihood of <3
            4. collapse overlaps 
            5. get the best breakpoints
        """
        print >>stderr, "resolving breakpoints..."
        final_calls = []
        
        for chr, overlapping_calls in self.overlapping_calls_by_chr.iteritems():
            print >>stderr, chr, "%d calls in this chr"%(len(overlapping_calls))
            
            indiv_cps = indiv_DTS.get_cps_by_chr(chr) 
            ref_cps = {}
            for ref, refDTS in ref_DTSs.iteritems():
                ref_cps[ref] = refDTS.get_cps_by_chr(chr) 
            
            g = genotyper(chr, gglob_dir=gglob_dir, plot_dir=out_viz_dir, subset_indivs = subset_indivs) 

            curr_dCGHs = self.get_curr_chr_dCGHs(chr, dCGHs)
            wnd_starts, wnd_ends = indiv_DTS.get_wnds_by_chr(chr)

            t=time.time()
            n_assessed=-1
            for overlap_cluster in overlapping_calls:
                overlap_cutoff = 0.85
                #overlap_cutoff = 0.75
                n_assessed+=1
                
                resolved_calls = overlap_cluster.overlap_resolve(overlap_cutoff, 
                                                                 ll_cutoff, 
                                                                 tbx_dups, 
                                                                 min_overlapping=min_overlapping) 
                
                if len(resolved_calls) == 0: continue
                
                variable_clusts = []
                for clust in resolved_calls:
                    """
                    now take all these resolved calls,  
                    and genotype to ensure this seg is var
                    """
                    d = self.get_delta(clust, wnd_starts, wnd_ends, curr_dCGHs, indiv_cps, ref_cps)
                    if clust.size == 1: 
                        print "skipping single call cluster"
                        continue
                    if d > 1.0:
                        variable_clusts.append(clust)
                    elif d>min_d:
                        X, idx_s, idx_e = g.get_gt_matrix(chr, clust.get_med_start(), clust.get_med_end())
                        gX = g.GMM_genotype(X)
                        if gX.gmm.n_components>1 and gX.is_var(indiv_id, g):
                            variable_clusts.append(clust)
                
                overlapping_call_clusts = get_overlapping_call_clusts(variable_clusts)

                for clust in overlapping_call_clusts:
                    final_call = self.get_final_call(clust)
                    final_calls.append(final_call)
                

                """
                if min(overlap_cluster.all_starts)>14454554 and min(overlap_cluster.all_starts)<14933135:
                self.plot_call(clust, wnd_starts, wnd_ends, curr_dCGHs, indiv_cps, ref_cps)
                raw_input()
                """
                #overlapping_call_clusts = get_overlapping_call_clusts(resolved_calls, flatten = True)
            print "n-final_calls:%d n_assessed:%d time_elapsed:%fs"%(len(final_calls), n_assessed, time.time()-t)

        print >>stderr, "done"
        print >>stderr, "%d calls with likelihood <%f"%(len(final_calls), ll_cutoff)
        return final_calls

    def get_final_call(self, resolved_calls):
        """
        the best call looks to be the min best start and max best end
        **COULD be an issue with DREADFUL long calls??? FOLLOW UP
        """
            
        mn = 1e100
        mx = -1

        ref_indivs = []
        contig = None

        for call in resolved_calls: 
            s,e = call.get_min_start_max_end()
            s = call.get_med_start()
            e = call.get_med_end()
            mn = min(s,mn)
            mx = max(e,mx)
            contig = call.chr
            ref_indivs+=call.ref_indivs
        
        ref_indivs=list(Set(ref_indivs))
        
        f_call = simple_call(contig, mn, mx, ref_indivs, resolved_calls)
        return f_call

    def get_delta(self, call_clust, wnd_starts, wnd_ends, dCGHs, indiv_cps, ref_cps):
        """
        the mean distance between an individual call and the refs it was called against
        """
        
        start, end = call_clust.get_med_start(), call_clust.get_med_end()
        w_s, w_e = np.searchsorted(wnd_starts, start), np.searchsorted(wnd_ends, end)
        
        ref_indivs = call_clust.ref_indivs
        ref_indivs=list(Set(ref_indivs))

        ds = []
        for ref in ref_indivs:
            d = np.mean(indiv_cps[w_s:w_e]-ref_cps[ref][w_s:w_e])
            ds.append(d)
        
        return np.mean(np.absolute(np.array(ds)))
        #print ds, np.mean(np.absolute(np.array(ds))), np.median(np.absolute(np.array(ds))) 
        

    def plot_call(self, resolved_calls, wnd_starts, wnd_ends, dCGHs, indiv_cps, ref_cps, fn_out="test.pdf"):
        
        """
        choose the call that maximizes the delta between the outside windows and the inside windows
        OR
        can i use the maximization algo below? delta = ((mu_middle-mu_left) + (mu_middle-mu_right))^2
        """
        plt.gcf().set_size_inches(14,6)
        fig, axes = plt.subplots(2,2) 
        
        font = {'family' : 'normal', 'weight': 'normal', 'size': 5}
        matplotlib.rc('font', **font) 
        
        colors = ['b','g','r','c','m','y','k']
        n_colors = len(colors)
        
        mn = 1e100
        mx = -1
        ref_indivs = []
        for call in resolved_calls: 
            s,e = call.get_min_start_max_end()
            print 'min,max', s,e
            print 'size,', call.size
            print "ll", call.get_log_likelihood()
            #s = call.get_best_start()
            #e = call.get_best_end()
            s = call.get_med_start()
            e = call.get_med_end()
            print 'best', s,e
            mn = min(s,mn)
            mx = max(e,mx)
            ref_indivs+=call.ref_indivs
        
        ref_indivs=list(Set(ref_indivs))
        w_s, w_e = np.searchsorted(wnd_starts,mn), np.searchsorted(wnd_ends,mx)
        
        #self.get_best_call(resolved_calls, wnd_starts, wnd_ends, ref_indivs, dCGHs)

        wnd_start = w_s - 20 
        wnd_end = w_e + 20 

        plt_s = wnd_starts[wnd_start]
        plt_e = wnd_ends[wnd_end]

        #t_vect = np.ones(wnd_end-wnd_start)
        t_vect = np.zeros(wnd_end-wnd_start)
        mids = (wnd_starts[wnd_start:wnd_end]+ wnd_ends[wnd_start:wnd_end])/2.0
        
        #t_vect = t_vect[:-2]
        #mids = mids[:-2]
         


        i =0
        for key, d in dCGHs.iteritems():
            #dsum =np.sum(d[wnd_start:wnd_end])
            #print dsum, colors[i%n_colors]
            #t_vect = t_vect*np.power(d[wnd_start:wnd_end],1)

            if key in ref_indivs:
                axes[0,0].plot(mids, d[wnd_start:wnd_end], colors[i%n_colors])
                cur_d = d[wnd_start:wnd_end]
                
                s,e, ps, pe = self.get_max_path(cur_d, wnd_starts[wnd_start:wnd_end], wnd_ends[wnd_start:wnd_end])
                mu = np.median(cur_d[ps:pe])
                print ps, pe
                print mu, s, e
                axes[0,0].plot([s,e], [mu,mu], colors[i%n_colors],marker='+', ls='-')
                axes[0,0].plot([s,e], [mu,mu], colors[i%n_colors])
                
                g2=ndi.gaussian_filter1d(cur_d,(wnd_end-wnd_start)/20.0,order=2)
                g1=ndi.gaussian_filter1d(cur_d,(wnd_end-wnd_start)/20.0,order=1)

                g2_10=ndi.gaussian_filter1d(cur_d,(wnd_end-wnd_start)/10.0,order=2)
                g1_10=ndi.gaussian_filter1d(cur_d,(wnd_end-wnd_start)/10.0,order=1)
                if np.sum(cur_d)<0:
                    cur_d = -1*cur_d
                
                t_vect = t_vect+cur_d

                #axes[0,1].plot(mids, g2, colors[i%n_colors])
                #axes[0,1].plot(mids, g2_10, colors[i%n_colors])
            else:
                axes[0,0].plot(mids, d[wnd_start:wnd_end], colors[i%n_colors], lw=.2)
            i+=1
        
        y=-.2
        for call in resolved_calls: 
            for indiv_call in call.calls:
                s,e = indiv_call['start'], indiv_call['end']
                axes[0,0].plot([s,e], [y,y], 'k', lw=1.5)
                y-=.1
        
        #t_vect = t_vect*t_vect
        axes[0,1].plot(mids, t_vect)
        axes[1,0].plot(mids, indiv_cps[wnd_start:wnd_end],'g', lw=1)
        axes[1,0].plot(mids, indiv_cps[wnd_start:wnd_end],'.g')
        for ref, cps in ref_cps.iteritems():
            axes[1,0].plot(mids, cps[wnd_start:wnd_end], lw=.1)
            
        for i in xrange(2):
            for j in xrange(2):
                axes[i,j].set_xlim([plt_s,plt_e])
                ym=axes[i,j].get_ylim()
                if i==0 and j ==1: ym = [1e-50,ym[1]]
                axes[i,j].plot([mn,mn],ym,'r')
                axes[i,j].plot([mx,mx],ym,'r')
                #print axes[i,j].get_yticklabels()
                #[i.set_fontsize(10) for i in axes[i,j].get_yticklabels()]

        axes[0,0].set_xlabel("indiv dCGH")
        axes[1,0].set_xlabel("cp")
        axes[0,1].set_xlabel("transformed_data")
        #axes[0,1].set_yscale("log")
        axes[1,1].set_xlabel("cp")
        plt.savefig(fn_out)
        plt.gcf().clear()
        
    
    def get_best_call(self, resolved_calls, wnd_starts, wnd_ends, ref_indivs, dCGHs):
        all_starts, all_ends = [], []
        for call in resolved_calls:
            all_starts += call.all_starts
            all_ends += call.all_ends
            
        all_starts = np.unique(np.array(all_starts)) 
        all_ends = np.unique(np.array(all_ends))

        min_start = all_starts[0]
        max_end = all_ends[-1]
        
        all_wss, all_wes = np.searchsorted(wnd_starts, all_starts), np.searchsorted(wnd_ends, all_ends)
        
        min_w_start = all_wss[0]
        max_w_end = all_wes[-1]

        if len(ref_indivs) < len(dCGHs.keys()): 
            """
            if there ARE controls
                get the mean of the controls... or median?
            """
            l  = max_w_end - min_w_start  
            n_controls = len(dCGHs)-len(ref_indivs)
            n_called = len(ref_indivs)
            a_controls = np.zeros((n_controls,l))
            a_called = np.zeros((n_called,l))

            i_controls, i_called = 0, 0
            
            for indiv, dCGH in dCGHs.iteritems(): 
                if indiv in ref_indivs:
                    a_called[i_called,:] = dCGH[min_w_start:max_w_end]
                    i_called+=1
                else:
                    a_controls[i_controls,:] = dCGH[min_w_start:max_w_end]
                    i_controls+=1
            
            med_controls = np.median(a_controls, 1)
            

            print med_controls, med_controls.shape
            print a_called
            return
        else:
            """
            if no controls
            """
            print "NO CONTROLS"
            print len(ref_indivs), len(dCGHs.keys())
            exit(1)
        
        #all_starts/all_ends sorted
        print all_starts
        print all_ends
        w_s, w_e = np.searchsorted(wnd_starts,mn), np.searchsorted(wnd_ends,mx)
        return


    def get_max_path(self, d, starts, ends):
        l = d.shape[0]
        return 0, l, 0, l
        if np.sum(d) <0:
            d = d*-1
        
        csum = np.cumsum(d)
        
        sums=np.ones((l,l))*1e10
        tb = np.zeros((l,l))

        """
            mu from i - j
            front seg and back seg
        """
        
        for j in xrange(10, l):
            sum_back_seg = csum[-1] - csum[j] 
            len_back_seg = l-j  
            for i in xrange(10,l):
                sum_front_seg = csum[i]
                len_front_seg = i
                mu = (sum_front_seg+sum_back_seg)/(len_back_seg+len_front_seg)
                sums[i,j] = np.sum(d[i:j]-mu)

        m = np.unravel_index(np.argmax(sums),sums.shape)
        return starts[m[0]], ends[m[1]], m[0], m[1]

    def filter_call_by_size_cp(self, final_call, indiv_cps, wnd_starts, wnd_ends, min_len=100000, cp_min=1.7, cp_max=2.3):
        """
        for calls > min_len
            if they look normal, kill them
        """
        s, e = final_call.start, final_call.end

        if e-s < min_len:
            return False
        
        wnd_start, wnd_end = np.searchsorted(wnd_starts,s), np.searchsorted(wnd_ends,e)
        med_cp = np.median(indiv_cps[wnd_start:wnd_end])

        if med_cp<cp_max and med_cp >cp_min:
            #if call is in the normal range (c < max c > min) then filter! 
            return True

        return False 

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
        
        curr_cluster = indiv_call_cluster(curr_chr) 

        overlapping_calls_by_chr[curr_chr].append(curr_cluster)
        n=1

        for idx, row in self.callset_table.pd_table.iterrows():
            chr, start, end = row['chr'], row['start'], row['end']
            if start<curr_max and curr_chr == chr:
                curr_cluster.add(row)
            else:
                curr_cluster = indiv_call_cluster(chr) 
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

def get_overlapping_call_clusts(calls):
    
    G = nx.Graph()
    
    for call in calls:
        G.add_node(call)

    for clust1 in calls:
        si, ei = clust1.get_med_start(), clust1.get_med_end()
        for clust2 in calls:
            sj, ej = clust2.get_med_start(), clust2.get_med_end()
            if intersect(si,ei,sj,ej):
                G.add_edge(clust1, clust2)
    return nx.connected_components(G)

class indiv_call_cluster:
    def __init__(self, chr):
        self.calls = []
        self.all_starts = []
        self.all_ends = []
        self.chr = chr
        self.size = 0
        self.log_likelihood =  None
        self.ref_indivs = []

    def get_min_start_max_end(self):
        min_s = min(self.all_starts)
        max_e = max(self.all_ends)
        return min_s, max_e

    def get_dup_overlap(self, tbx_dups):
        #assuming dup file is collapsed
        min_s = min(self.all_starts)
        max_e = max(self.all_ends)
        t = 0
            
        for l in tbx_dups.fetch(self.chr,min_s,max_e,parser=pysam.asTuple()):
            c,s,e = l
            s,e = int(s), int(e)
            curr_s = s>min_s and s or min_s
            curr_e = e<max_e and e or max_e
            t += curr_e-curr_s
        if t==0: return 0.0
        
        return float(t) / float(max_e - min_s)

    def add(self, call):
        self.calls.append(call)
        self.all_starts.append(call['start'])
        self.all_ends.append(call['end'])
        self.ref_indivs.append(str(call['indiv_ref']))
        self.size+=1
    
    def over_simplify(self):
        """
        take all calls w/ the same start, and average their ends
        """
        print >>stderr, "over-simplifying..."
        curr_call=[]
        new_calls = []
        new_all_starts = []
        new_all_ends = []
        curr_ends = [] 
        for call in self.calls:
            if len(curr_call)==0:
                curr_call=call
                new_calls.append(curr_call)
                new_all_starts.append(curr_call['start'])
                curr_ends = [curr_call['end']]
                continue
            if call['start'] == curr_call['start']:
                curr_call['p'] = 10**(np.log(curr_call['p'])+np.log(call['p'])) 
                curr_ends.append(curr_call['end'])
            else:
                end = int(np.median(np.array(curr_ends)))
                curr_call['end'] = end
                new_all_ends.append(curr_call['end'])
                
                curr_call=call
                new_calls.append(curr_call)
                new_all_starts.append(curr_call['start'])
                curr_ends = [curr_call['end']]
        
        self.calls=new_calls
        self.size=len(new_calls)
        self.all_starts = new_all_starts
        self.all_ends = new_all_ends


    def simplify(self):
        """
        take calls w/ identical breakpoints and pre-merge them 
        reduces the complexity of the next n^2 operation
        """
        curr_call=[]
        new_calls = []
        new_all_starts = []
        new_all_ends = []
        
        for call in self.calls:
            if len(curr_call)==0:
                curr_call=call
                new_calls.append(curr_call)
                new_all_starts.append(curr_call['start'])
                new_all_ends.append(curr_call['end'])
                continue
            if call['start'] == curr_call['start'] and call['end'] == curr_call['end']:
                curr_call['p'] = 10**(np.log(curr_call['p'])+np.log(call['p'])) 
            else:
                curr_call=call
                new_calls.append(curr_call)
                new_all_starts.append(curr_call['start'])
                new_all_ends.append(curr_call['end'])
        
        self.calls=new_calls
        self.size=len(new_calls)
        self.all_starts = new_all_starts
        self.all_ends = new_all_ends


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
    
    def get_call_str(self):
        ret_str = ""
        call = self.calls[0]
        
        size = np.median(np.array([ c['window_size'] for c in self.calls]))
        mu = np.median(np.array([ c['mu'] for c in self.calls]))
        
        ret_str="\t".join([ str(s) for s in  [ call["indiv_ref"], 
                                               call["indiv_test"], 
                                               call["chr"],
                                               self.get_best_start(), 
                                               self.get_best_end(), 
                                               mu,
                                               10.0**self.get_log_likelihood(),
                                               size ]])
        return ret_str

    def get_bed_str(self):
        
        call = self.calls[0]
        ret_str="%s\t%d\t%d\t%f\t%f\t+\t%d\t%d\t0,0,0"%(call["chr"],
                                                     self.get_best_start(),
                                                     self.get_best_end(),
                                                     10.0**self.get_log_likelihood(),
                                                     10.0**self.get_log_likelihood(),
                                                     self.get_best_start(),
                                                     self.get_best_end())
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

    def get_med_start(self):
        return np.median(np.array(self.all_starts))
    
    def get_med_end(self):
        return np.median(np.array(self.all_ends))

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

        plt.savefig('test2.png')
        plt.gcf().clear()
        

        
    def cluster_by_recip_overlap(self, cutoff, min_overlapping, ll_cutoff=0, plot=False):
        """
        takes the current cluster - a set of calls each w/ llls, and 
        clusters them based on their recdip overlap
        returns a new set of clusters made up of original set of clusters
            - make sure they pass ll cutoff
            - make sure pass min_size cutoff
        """

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
        
        new_groups = []
        for grp in np.unique(grps):
            new_group = indiv_call_cluster(self.chr) 
            for idx in np.where(grps==grp)[0]:
                new_group.add(self.calls[idx])
            if new_group.get_log_likelihood()<ll_cutoff and new_group.size >= min_overlapping: 
                new_groups.append(new_group)
        return new_groups
        
    def overlap_resolve(self, overlap_cutoff, ll_cutoff, tbx_dups, min_overlapping=1, do_plot=False):
        """
        though a bunch of calls may overlap, they may not be all 'good'
        or all the same call. Thus, further cluster overlapping calls
        into those with reciprocal overlap properties
            
        ***testing, output just recip overlap calls
        criteria for a call
        1. at least 2 calls reciprical overlap and look the same
        2. a.) after recip overlpa the sum of the lls is <ll_cutoff (-ve)
        2. b.) after recip overlap still at least 2 calls
        """
        
        frac_dup = self.get_dup_overlap(tbx_dups)
        
        #recall, SIZE means the number of events in the cluster
        if self.size <min_overlapping: 
            return []

        #if self.size>1000: 
        #    self.over_simplify()
        #else:
        #    self.simplify()
        
        """
        single event (in the event that we ARE allowing them, by default)
        OR, if by simplify is now just one call
        """
        if self.size == 1:
            #return self.get_log_likelihood() < ll_cutoff and [CNV_call([self],self.chr)] or []
            #return self.get_log_likelihood() < ll_cutoff and [self] or []
            recip_overlap_clusts = [self]
        else: 
            recip_overlap_clusts = self.cluster_by_recip_overlap(overlap_cutoff, 
                                                                 min_overlapping=min_overlapping, 
                                                                 ll_cutoff=ll_cutoff, 
                                                                 plot=do_plot)
        return recip_overlap_clusts



class call_cluster(object):
    """
    similar to indiv_call_cluster but for use on
    calls made amongst individuals
    """
    
    def __init__(self):
        self.calls = []
        self.starts = []
        self.ends = []
        self.lengths = []
        
    def add(self, call):
        self.calls.append(call)
        self.starts.append(call['start'])
        self.ends.append(call['end'])
        self.lengths.append(call['end']-call['start'])
    
    def get_min_max_start_end(self):
        return min(self.starts), max(self.ends)
    
    def get_med_start_end(self):
        l = len(self.calls)
        idx = l/2
        return sorted(self.starts)[idx], sorted(self.ends)[idx]



    def frac_overlap(self, i, j, max_dif=-1, max_frac_uniq=-1):
        si, ei = self.starts[i],  self.ends[i]
        sj, ej = self.starts[j],  self.ends[j]

        if not (si < ej and sj < ei): 
            return 0   

        o = min(ei,ej) - max(si, sj)
        
        if max_dif!=-1 and max((ei-si)-o, (ej-sj)-o) > max_dif:
            return 0
        
        li=float(ei-si)
        lj=float(ej-sj)
        
        if max_frac_uniq!=-1 and max(1-(o/li), 1-(o/lj)) > max_frac_uniq:
            return 0

        return min(float(o)/(ei-si), 
                   float(o)/(ej-sj))
    
    def cluster_by_recip_overlap(self, cutoff, max_dif=-1, max_frac_uniq=-1):
        l = len(self.calls)
        mat = np.zeros((l,l))

        for i in xrange(l):
            for j in xrange(l):
                mat[i,j] = 1-(self.frac_overlap(i,j, max_dif, max_frac_uniq))
        
        mat=1.-dist.squareform( (1-mat) * (np.diag(np.repeat(1,l),0) -1) * -1)
        Z = hclust.linkage(mat, method='complete')
        grps = hclust.fcluster(Z, cutoff, criterion='distance')
        
        new_groups = []
        for grp in np.unique(grps):
            new_group = call_cluster()
            for idx in np.where(grps==grp)[0]:
                new_group.add(self.calls[idx])
            new_groups.append(new_group)

        return new_groups
        
def linkage_cluster(mat, cutoff):
    
    l = mat.shape[0]
    mat=1.-dist.squareform( (1-mat) * (np.diag(np.repeat(1,l),0) -1) * -1)
    Z = hclust.linkage(mat, method='complete')
    grps = hclust.fcluster(Z, cutoff, criterion='distance')
    return grps 


def seg_sets_intersect(segs1, segs2, max_frac_uniq = -1):
    """
    two sets of segs in an individual, we want to cluster
    come up w/ some distance
    """
    t1 = np.sum(np.array([s_e[1]-s_e[0] for s_e in segs1]))
    t2 = np.sum(np.array([s_e[1]-s_e[0] for s_e in segs2]))
    
    o=0 
    for se1 in segs1:
        for se2 in segs2:
            s1, e1 = se1
            s2, e2 = se2
            if s1<=e2 and s2<=e1:
                o += min(e1,e2)-max(s1,s2)
    
    o=float(o)
    min_frac_o = min(o/t1, o/t2)
    if max_frac_uniq!=-1 and 1-min_frac_o>=max_frac_uniq:
        return 0
    return min_frac_o


class cluster_callsets(object):
    """
    a class to cluster calls from multiple individuals
    into sets fot genotyping attempting to separate out
    calls that are overlapping but very different in their 
    breakpoints into distinct clusters
    """
    
    @classmethod
    def plot(cls, overlapping_call_clusts, out_dir, g, indivs_by_seg, s_es, CNV_s_e, cnv_segs_by_indiv):


        """
        row col
        0, 0 are the individual calls, clustered together
        1, 0 are the median s e of the calls AND the chopped up bits of those
        0, 1 is takign all the the CNV segs and clustering them by indiv
        1, 1 is taking those clustered individuals and the outputting the blocks
        """
         
        min_start = 9e9
        max_end = -1
        j=0
        for call_clust in overlapping_call_clusts:
            #print "recip overlap cluster %d, contains %d calls"%(j, len(call_clust.starts))
            min_start = min(min_start, min(call_clust.starts))
            max_end = max(max_end, max(call_clust.ends))
            j+=1

        print min_start, max_end, len(overlapping_call_clusts)
            
        plt.gcf().set_size_inches(14,6)
        fig, axes = plt.subplots(2,2) 
        
        font = {'family' : 'normal', 'weight': 'normal', 'size': 5}
        matplotlib.rc('font', **font) 
        
        colors = ['b','g','r','c','m','y','k']
        n_colors = len(colors)
        
        y=-.2
        y2=-.2
        #the individual clusters
        for j, call_clust in enumerate(overlapping_call_clusts):

            axes[1,0].plot([np.median(call_clust.starts), np.median(call_clust.ends)], [y2,y2], colors[j%n_colors], lw=1.5)
            y2-=.2
            for i, start in enumerate(call_clust.starts):
                end =  call_clust.ends[i]
                axes[0,0].plot([start, end], [y,y], colors[j%n_colors], lw=1.5)
                y-=.1
            y-=.2
        

        #all the little sub-bits
        y2-=.5
        for s_e in s_es:
            axes[1,0].plot(s_e, [y2,y2], 'r', lw=2.5)
            y2-=.2

        y2-=.2
        for s_e in CNV_s_e:
            axes[1,0].plot(s_e, [y2,y2], 'g', lw=2.5)
            y2-=.1
        
        ###############
        
        i=0
        y3=0
        for seg, indivs in indivs_by_seg.iteritems():
            for indiv in indivs:      
                axes[1,1].plot(seg,[y3,y3], colors[i%n_colors], lw=1.5)
                y3-=1
            i+=1
        #p = axes[1,1].pcolor(c)
        #fig.colorbar(p, cax=axes[0,1])

        y3=0
        i=0 
        #each indivs cnv segs
        for indiv, cnv_segs in cnv_segs_by_indiv.iteritems():
            for seg in cnv_segs:
                axes[0,1].plot(seg,[y3,y3], colors[i%n_colors], lw=1.5)
            y3-=.1
            i+=1

        axes[0,0].set_xlabel("calls")
        axes[0,0].set_xlim([min_start, max_end])
        axes[0,0].set_ylim([0, y-.2])
        
        axes[1,0].set_xlim([min_start, max_end])
        axes[1,0].set_ylim([0, y2-.2])
        axes[0,1].set_xlim([min_start, max_end])
        
        ylims = axes[1,1].get_ylim()
        axes[1,1].set_ylim([ylims[0]-1, ylims[1]+1])

        axes[1,1].set_xlim([min_start, max_end])
        plt.savefig("%s/%d_%d_bps.pdf"%(out_dir, min_start, max_end))
        plt.gcf().clear()
        

    def __init__(self, fn_table, contig):
        
        self.contig = contig
        self.callset_table = simple_callset_table(fn_table)
        self.callset_table.sort()
        self.callset_table.filter_by_chr(self.contig)
    
    #def get_genotypable_loci(self, total_subsets, subset):

    def get_overlapping_call_clusts(self, total_subsets, subset):
        """
        get sets of overlapping calls and return a list
        of those calls grouped into recip-overlapping clusters
        """
        lists_of_overlapping_calls, n = self.get_overlapping_calls()
        print >>stderr, "%d overlapping calls"%n
        
        for overlapping_call_cluster in self.subset_call_list(lists_of_overlapping_calls, total_subsets, subset):
            if len(overlapping_call_cluster.calls)==1:
                yield [overlapping_call_cluster]
            else:
                yield overlapping_call_cluster.cluster_by_recip_overlap(0.35, max_dif=15000)
                #yield overlapping_call_cluster.cluster_by_recip_overlap(0.35)

    def subset_call_list(self, overlap_call_groups, total_subsets, subset):
        
        l = len(overlap_call_groups)
        call_groups_to_assess = []

        i=subset
        while i<l:
            call_groups_to_assess.append(overlap_call_groups[i])
            i+=total_subsets
        
        print >>stderr, "%d calls in subset"%(len(call_groups_to_assess))
        return call_groups_to_assess


    def get_overlapping_calls(self):
        """
        glob up all those calls that physically overlap
        assumes self.callset_table is sorted
        """
        curr_max=self.callset_table.pd_table.iloc[0]['end']
        
        overlapping_calls =  []
        curr_cluster = call_cluster() 

        overlapping_calls.append(curr_cluster)
        n=1

        for idx, call in self.callset_table.pd_table.iterrows():
            chr, start, end = call['chr'], call['start'], call['end']
            if start<curr_max:
                curr_cluster.add(call)
            else:
                curr_cluster = call_cluster() 
                curr_max = end
                curr_cluster.add(call)
                overlapping_calls.append(curr_cluster)
                n+=1
            
            curr_max = max(end, curr_max) 
        
        return overlapping_calls, n
        


