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

from sets import Set

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

    def filter(self, p_cutoff, size_cutoff, single_window_cutoff, divide_by_mu=False):
        """
        apply size cutoff (no calls of size 1 window!)
        """
        print >>stderr, "parsing calls with p<=%f and l>=%d..."%(p_cutoff, size_cutoff)
        if divide_by_mu: 
            self.pd_table = self.pd_table[((self.pd_table['p']/np.exp((self.pd_table['mu']**2)))<p_cutoff)]
        else:
            self.pd_table = self.pd_table[(self.pd_table['p']<p_cutoff)]

        self.pd_table = self.pd_table[(self.pd_table['window_size']>=size_cutoff)]
        ##for windows of size 1, if they exist, ensure they exceed the max_mu
        #self.pd_table = self.pd_table[((self.pd_table['window_size']>1)|(np.absolute(self.pd_table['mu'])>=single_window_cutoff))]
        #self.pd_table = self.pd_table[( (self.pd_table['window_size']>=size_cutoff) |
        #                                (self.pd_table['p']>0) )]
        print >>stderr, "done"


class simple_call:
        
    def __init__(self, contig, start, end, ref_indivs, resolved_calls):
        self.contig = contig
        self.start = start
        self.end = end
        self.ref_indivs  = ref_indivs #indivs CALLED over region
        self.resolved_calls = resolved_calls

        self.total_ll = np.sum(np.array([call.get_log_likelihood() for call in self.resolved_calls]))


class indiv_cluster_calls:
    
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


    def resolve_overlapping_clusters(self, ll_cutoff, tbx_dups, indiv_DTS, dCGHs, verbose=False, min_overlapping=2):
        """
        resolve overlapping clusters into individual calls 
        let the calls have likelihoods
            1. do the recip overlap cluster - clusters very similar calls w/ similar break-points
            2. make sure there are at least 2 calls in there (w/ similar breakpoints)
            *3. also make sure those calls are reasonable (not cp 2) ???really????
            4. collapse overlaps 
            5. get the best breakpoints
        """
        print >>stderr, "resolving breakpoints..."
        final_calls = []
        
        for chr, overlapping_calls in self.overlapping_calls_by_chr.iteritems():
            print >>stderr, chr
            
            indiv_cps = indiv_DTS.get_cps_by_chr(chr) 
            curr_dCGHs = self.get_curr_chr_dCGHs(chr, dCGHs)
            wnd_starts, wnd_ends = indiv_DTS.get_wnds_by_chr(chr)

            for overlap_cluster in overlapping_calls:
                overlap_cutoff = 0.85
                
                resolved_calls = overlap_cluster.overlap_resolve(overlap_cutoff, 
                                                                 ll_cutoff, 
                                                                 tbx_dups, 
                                                                 min_overlapping=min_overlapping) 
                
                if len(resolved_calls) == 0: continue
                #overlapping_call_clusts = get_overlapping_call_clusts(resolved_calls, flatten = True)

                final_call = self.get_final_call(resolved_calls)
                """
                self.plot_call(resolved_calls, wnd_starts, wnd_ends, curr_dCGHs, indiv_cps)
                raw_input()
                """
                #overlapping_call_clusts = get_overlapping_call_clusts(resolved_calls, flatten = True)
                if not self.filter_call_by_size_cp(final_call, indiv_cps, wnd_starts, wnd_ends):
                    final_calls += [final_call]

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
            s = call.get_best_start()
            e = call.get_best_end()
            mn = min(s,mn)
            mx = max(e,mx)
            contig = call.chr
            ref_indivs+=call.ref_indivs
        
        ref_indivs=list(Set(ref_indivs))
        
        f_call = simple_call(contig, mn, mx, ref_indivs, resolved_calls)
        return f_call

    def plot_call(self, resolved_calls, wnd_starts, wnd_ends, dCGHs, indiv_cps):
        
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
            s = call.get_best_start()
            e = call.get_best_end()
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
        axes[1,0].plot(mids, indiv_cps[wnd_start:wnd_end])
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
        plt.savefig('test.pdf')
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

def get_overlapping_call_clusts(calls, flatten = False):
    
    G = nx.Graph()
    
    for call in calls:
        G.add_node(call)

    for clust1 in recip_overlap_clusts:
        si, ei = clust1.get_best_start(), clust1.get_best_end()
        for clust2 in recip_overlap_clusts:
            sj, ej = clust2.get_best_start(), clust2.get_best_end()
            if intersect(si,ei,sj,ej):
                G.add_edge(clust1, clust2)
    if not flatten:    
        return nx.connected_components(G)
    else:
        flat_clusts = []
        for overlap_c in nx.connected_components(G):
            start = overlap_c[0].get_best_start()
            end = overlap_c[0].get_best_end()
            ll = 0 
            for c in overlap_c:
                if c.get_best_end() >end: 
                    end = c.get_best_end()
                if c.get_best_start() <start: 
                    end = c.get_best_end()
                ll+=c.get_log_likelihood()
            curr_clust = call_cluster(overlap_c[0].chr)
             
            curr_clust.add({"chr":curr_clust.chr,
                            "start":start,
                            "end":end,
                            "p":10**ll,
                            "window_size":-1,
                            "mu":-1,
                            "indiv_ref":"flattened",
                            "indiv_test":"flattened"})
            flat_clusts.append([curr_clust])

        return flat_clusts

class call_cluster:
    def __init__(self, chr):
        self.calls = []
        self.all_starts = []
        self.all_ends = []
        self.chr = chr
        self.size = 0
        self.log_likelihood =  None
        self.ref_indivs = []

        """
        'complex' refers to calls that have been 
        clustered (overlap) together, but reach significance
        self.compound = False 
        """
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
        self.ref_indivs.append(call['indiv_ref'])
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
            raw_input()
        
        new_groups = []
        for grp in np.unique(grps):
            new_group = call_cluster(self.chr) 
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

        criteria for a call
        1. at least 2 calls reciprical overlap and look the same
        2. a.) after recip overlpa the sum of the lls is <ll_cutoff (-ve)
        2. b.) after recip overlap still at least 2 calls
        """
        
        frac_dup = self.get_dup_overlap(tbx_dups)
        
        #recall, SIZE means the number of events in the cluster
        if self.size <min_overlapping: 
            return []

        if self.size>1000: 
            self.over_simplify()
        else:
            self.simplify()
        
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


class cluster_
