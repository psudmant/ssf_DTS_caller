from traverse_contours import get_contours
#from edge_cluster import *
from c_hierarchical_edge_merge import *

import numpy as np
from sys import stderr
import scipy.ndimage as ndi

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

import time
import pysam
import bisect

from null_distribution import null_distribution
from get_windowed_variance import get_windowed_variance

class call:
    def __init__(self, chr, start, end, value, wnd_start, wnd_end, values):

        self.chr = chr
        self.start = start
        self.end = end
        self.value = value
        self.wnd_start = wnd_start
        self.wnd_end = wnd_end
        
        self.width = self.wnd_end - self.wnd_start + 1
         
        self.values = values
        
        self.p_value = None
        self.fdr_significant = False
        self.significance_level = -1
        self.t_stats = None


class callset:
    
    def __init__(self):
        self.calls = []
    
    def add_call(self, chr, start, end, val, wnd_end, wnd_start, values):
        self.calls.append(call(chr, start, end, val, wnd_end, wnd_start, values))

    def __iadd__(self, other):  
        self.calls += other.calls
        return self
    
    def get_t_stats(self,null_dist):
        print "getting t-values..." 
        t = time.time() 
        for call in self.calls:
            call.t_stats = null_dist.get_t_stats(call) 
        print "done in %fs"%(time.time() - t)
        
    def get_p_values(self,null_dist):

        print "getting p-values..." 
        t = time.time() 
        for call in self.calls:
            #call.p_value = null_dist.get_t_test_p_value(call) 
            call.p_value = null_dist.get_corrected_t_test_p_value(call) 

            #call.p_value = null_dist.get_rank_sum_p_value(call) 
            #call.p_value = null_dist.get_mu_based_p_value(call) 
            #call.p_value = null_dist.get_ll_based_p_value(call) 
        print "done in %fs"%(time.time() - t)
            
    def get_bh_corrected_significant(self,fdr):
        """
        Benjamini Hochberg correction for 
        significance at some fdr
        """
        print "correcting p-values by Benjamini Hochberg procedure at fdr %f..."%fdr
        self.calls = sorted(self.calls,key=lambda x: x.p_value)
        
        m = float(len(self.calls))
        
        max_k = 0

        for k, call in enumerate(self.calls):
            call.significance_level = k 
            if call.p_value < (fdr * (k+1)/m):
                max_k = k

        for k in xrange(max_k+1):
            self.calls[k].fdr_significant = True
        
        print "done"
    
    def get_calls_in_range(self,wnd_start,wnd_end):
        calls = np.array([ call for call in self.calls if (call.wnd_end>wnd_start and call.wnd_start<wnd_end ) ])
        return calls
        
    def output(self,fn_out,t_stats=False):
        
        FOUT = open(fn_out,'w')
        additional=""
        for i,call in enumerate(sorted(self.calls,key=lambda x:x.p_value)):
            if t_stats:
                additional = "\t%(mu1)f\t%(var1)f\t%(n1)d\t%(null_mu)f\t%(null_var)f\t%(null_n)d"%(call.t_stats)

            print >>FOUT, "%d\t%s\t%d\t%d\t%f\t%.8e\t%s\t%d%s"%(i,
                                                    call.chr, 
                                                    call.start, 
                                                    call.end,
                                                    call.value,
                                                    call.p_value,
                                                    call.fdr_significant and "*" or "-",
                                                    (call.wnd_end-call.wnd_start),
                                                    additional)

class ssf_caller:
    def __init__(self,chr,cp_data, starts, ends, cutoff_scale, **kwargs):
        
        max_merge=kwargs.get("max_merge",0.5)
        use_means=kwargs.get("use_means",False)
        n_scales=kwargs.get("n_scales",51)
        #n_scales=kwargs.get("n_scales",30)
        scale_width=kwargs.get('scale_width',1)
        n_bin_smoothings=kwargs.get('-n_bin_smoothings',0)
        smoothing_kernel=kwargs.get('smoothing_kernel',np.array([1,2,1]))
        self.chr = chr
        self.cutoff_scale =  cutoff_scale
        self.scales = list(np.arange(1,n_scales,scale_width))
        self.starts = starts
        self.ends = ends
        self.n_wnds = self.starts.shape[0]
        self.cp_data = cp_data
        
        self.der1=np.zeros((len(self.scales),self.n_wnds),dtype=np.float32)
        self.der2=np.zeros((len(self.scales),self.n_wnds),dtype=np.int8)
        
        self.vars = get_windowed_variance(cp_data.astype(np.float64),500) 
        self.l_vars = np.roll(self.vars,501) 
        self.r_vars = np.roll(self.vars,-501) 

        print >>stderr, "scales range from %f-%f"%(self.scales[0],self.scales[-1])
        for i in xrange(n_bin_smoothings):
            print >>stderr,"doing binomial smooth #%d"%i
            cp_data=ndi.convolve1d(cp_data,smoothing_kernel)/np.sum(smoothing_kernel)
        
        transitions_by_scale = {}
        print >>stderr, "finding contours..."
        for i_scale,scale in enumerate(self.scales):
            stderr.write("%.2f "%(scale))
            stderr.flush()
            g1=ndi.gaussian_filter1d(cp_data,scale,order=1)                            
            g2=ndi.gaussian_filter1d(cp_data,scale,order=2)
            edges,pos_edges,neg_edges = self.get_n_edges(g1,g2)
            self.der1[i_scale,:]=g1
            self.der2[i_scale,:]=pos_edges-neg_edges
            transitions_by_scale[scale]=(edges,pos_edges,neg_edges)
        stderr.write("done\n")
        
        self.contour_intersects,x_intercept_to_scale=get_contours(self.der2)    
        
        ######NOW we have all the per-scale contours
        #print contour_intersects
        edges_passing_cutoff =[]
        curr_all_edges=[]
        curr_all_edges_scales=[]
       
        #take all the edges discovered at some scale
        for scale,edges in self.contour_intersects.iteritems():
            curr_all_edges.extend(edges)
            curr_all_edges_scales.extend([scale for i in xrange(len(edges))])
            if scale >=cutoff_scale:
                edges_passing_cutoff.extend(edges)
        edges_passing_cutoff=sorted(set(edges_passing_cutoff))
            
        all_edges_scales=sorted(zip(curr_all_edges,curr_all_edges_scales))
        stderr.write("hierarchically merging segments\n")
        
        t = time.time()
        segments_s, segments_e, cps = c_hierarch_merge_edges(cp_data, 
                                                        edges_passing_cutoff,
                                                        max_merge,
                                                        use_means,
                                                        self.n_wnds,
                                                        self.starts,
                                                        self.ends)
        #segments_s, segments_e, cps = hierarch_merge_edges(cp_data, 
        #                                                edges_passing_cutoff, 
        #                                                max_merge,use_means)
        self.segment_edges=(segments_s,segments_e,cps)
        print >>stderr, "hierarchical clustering completed in %fs"%(time.time()-t)  
    
    def get_exclude_coords(self, ex_starts, ex_ends):
                                                    
        mx=self.starts.shape[0]-1
        n_exclude = len(ex_ends)     
        ex_wnd_starts = np.searchsorted(self.starts, ex_starts)
        ex_wnd_ends   = np.searchsorted(self.ends, ex_ends)
        ex_wnd_starts = np.amax(np.c_[ex_wnd_starts-1,np.zeros(n_exclude)],1).astype(int)
        ex_wnd_ends = np.amin(np.c_[ex_wnd_ends+1,np.ones(n_exclude)*mx],1).astype(int)
        ex_starts = self.starts[ex_wnd_starts] 
        ex_ends = self.ends[ex_wnd_ends] 

        ex_coords = [] 
        
        curr_s = ex_starts[0]
        curr_e = ex_ends[0]

        #print ex_wnd_starts
        #print ex_wnd_ends

        for i in xrange(1, n_exclude):
            if ex_starts[i] < curr_e:
                curr_e = ex_ends[i]
            else:
                ex_coords.append(tuple([curr_s,curr_e]))
                curr_s = ex_starts[i]
                curr_e = ex_ends[i]
        
        ex_coords.append(tuple([curr_s,curr_e]))
        return ex_coords

    def subtract_excluded(self, wnd_start, wnd_end, ex_coords):
        """
        subtract the exclusion coordinates from the full call
        ea. output '.'s exclude 'x's over cal '-'s
           ------------------------
          xxxx....xx....xxx....xxxxxx
        """
        start = self.starts[wnd_start]
        end = self.ends[wnd_end]

        #wnd_start-wnd_end totally encompassed by gap
        if start >= ex_coords[0][0] and end <= ex_coords[0][1]:
            return []
        
        #starting wnd, either the end of the first gap, or the start of the wnd
        if start >= ex_coords[0][0]:
            init = ex_coords[0][1] 
            ex_coords.pop(0)
        else:
            init = start
        
        #ending wnd, either the start of the last gap, or the end of the wnd
        if len(ex_coords)>0 and end<=ex_coords[-1][1]:
            final = ex_coords[-1][0]  
            ex_coords.pop()
        else:
            final = end 
        
        regions = [init, final]
        
        for ex in ex_coords: 
            regions.append(ex[0])
            regions.append(ex[1])
        
        regions = sorted(regions)
        
        final_wnds = []
        
        i=0
        while i<len(regions):
            final_wnds.append(tuple([np.where(self.starts == regions[i])[0][0],
                                     np.where(self.ends == regions[i+1])[0][0]
                                    ]))
            i+=2
        
        return final_wnds
        

    def get_callset(self, exclude_tbxs=[], min_exclude_ratio=0.3, min_exclude_len=20000):
        """
        return segments and their copies in genome 
        coordinates
        adding subtraction of gaps
        """
        c=callset()

        wnd_starts,wnd_ends,cps = self.segment_edges
        wnd_starts,wnd_ends,cps = np.array(wnd_starts), np.array(wnd_ends), np.array(cps)
        
        for i in xrange(len(wnd_starts)-1):
            start, end = self.starts[wnd_starts[i]], self.ends[wnd_ends[i]]
            wnd_start, wnd_end = wnd_starts[i], wnd_ends[i]
            
            #exclude totally anything in these tbxs
            for exclude_tbx in exclude_tbxs:
                ex_starts, ex_ends = [], []
                for l in exclude_tbx.fetch(self.chr,start,end,parser=pysam.asTuple()):
                    _chr,_s,_e = l
                    _s, _e = int(_s), int(_e)
                    if _e-_s > min_exclude_len:
                        ex_starts.append(_s)
                        ex_ends.append(_e)

            n_exclude = len(ex_starts) 
            if n_exclude:
                ex_coords = self.get_exclude_coords(ex_starts, ex_ends)

                wnd_start_ends = self.subtract_excluded(wnd_start, wnd_end, ex_coords)
            else:
                wnd_start_ends = [tuple([wnd_start, wnd_end])]
            
            for i in xrange(len(wnd_start_ends)):
                wnd_start = wnd_start_ends[i][0]
                wnd_end = wnd_start_ends[i][1]

                c.add_call(self.chr,
                           self.starts[wnd_start],
                           self.ends[wnd_end],
                           np.mean(self.cp_data[wnd_start:wnd_end]),
                           wnd_start,
                           wnd_end,
                           self.cp_data[wnd_start:wnd_end])
        return c 

    def get_n_edges(self,der1,der2,n=0):
        """
        get the edges from the derivatives
        der1 is magnitudes        
        der2 is 0 crossings
        d_der2=diff(der2)
                 _________________________
                /                         \
               /                           \
         ------                             --------
        d1   0 to +ve                    0 to -ve 
        d2   +ve to 0 -ve                -ve to 0 +ve
        pm   -ve                           +ve
        ff   -ve (and should be +1)        +ve 
        
        pos to a neg (left to right) results in a -ve pm
        neg to a pos (left to right) results ina  +ve pm
         +1+1+1-1-1-1+1+1+1
         0 0-2 0 0+2 0 0
        """
        
        if n==0:
            n=der1.shape[0]    
        der2_pos=der2>0
        der2_neg=der2<0
        
        pm_array=np.zeros(der2.shape[0])
        pm_array=pm_array+der2_pos-der2_neg
        
        diffarray=np.diff(pm_array)
        
        #GET INNER THEN OUTER
        """
        this means the first bp is INSIDE the variant
        on the left hand side,
        on the right side the call is outside the variant 
        """
        intercepts_rise=np.where(diffarray<0)[0]+1
        intercepts_fall=np.where(diffarray>0)[0]+1
        #GET INNER COORDS
        #intercepts_rise=np.where(diffarray<0)[0]+1
        #intercepts_fall=np.where(diffarray>0)[0]
        #GET OUTER COORDS
        #intercepts_rise=np.where(diffarray<0)[0]
        #intercepts_fall=np.where(diffarray>0)[0]+1
        
        all_intercepts = np.unique(np.r_[intercepts_rise, 
                                         intercepts_fall,
                                         0,pm_array.shape[0]-1])
        intercepts=all_intercepts
        
        #magnitudes at each intercept
        mags = np.abs(der1[intercepts]) 
        
        len_mags=mags.shape[0]
        #locations in the intercepts array of the highest mags
        top_sorted_inds=np.argsort(mags)[(len_mags-n):len_mags] 
        
        #locations in the intercepts array of the highest mags
        #top_sorted_inds=np.argsort(mags) 
        edge_array=np.zeros(der1.shape[0])

        edge_array[intercepts[top_sorted_inds]]=1
        neg_edges=edge_array*-1*der2_neg    
        pos_edges=edge_array*der2_pos
        return edge_array,neg_edges,pos_edges

