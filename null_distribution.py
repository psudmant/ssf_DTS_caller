import numpy as np
import time
import pysam
from get_windowed_variance import get_windowed_variance

import scipy.stats as stats

def plot2d(density,x_edges, y_edges):
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   density=np.log(density) 
   plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
   #fig = plt.figure()
   fig,ax = plt.subplots()
   m = np.amax(density)
   #fig.set_figwidth(11)
   #fig.set_figheight(8.5)
   #box = [.05,.05,.9,.9]
   #ax = fig.add_axes(box)
   #ax, c = plt
   #print density
   extent=[y_edges[0],y_edges[-1],x_edges[1],x_edges[-2]]
   cax=ax.imshow(density, extent=extent, aspect='auto',interpolation='nearest')
   #cax=ax.imshow(density, extent=extent, aspect='auto')
   fig.colorbar(cax,orientation='horizontal',ticks=[0,m])
   #ax.colorbar()
   plt.savefig("./density.pdf",format="pdf")

def plot_density(density, edges, f):
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   
   plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
   fig = plt.figure()
   fig.set_figwidth(11)
   fig.set_figheight(8.5)
   box = [.05,.05,.9,.9]
   ax = fig.add_axes(box)
   ax.plot(edges[:-1],density,color='yellow')
   ax.fill_between(edges[:-1],0,density,alpha=.2,color='red')
   #ax.set_yscale('symlog')
   ax.set_xlim(edges[1],edges[-2])
   fig.savefig(f,format="pdf")

def plot_cum_density(density, edges, where):
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   
   plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
   fig = plt.figure()
   fig.set_figwidth(11)
   fig.set_figheight(8.5)
   box = [.05,.05,.9,.9]
   ax = fig.add_axes(box)
   ax.plot(edges[:-1],density,color='yellow')
   ax.fill_between(edges[:-1],0,density,alpha=.2,color='red')
   
   ax.fill_between(edges[:where],0,density[:where],alpha=.2,color='blue')

   ax.set_xlim(edges[1],edges[-2])
   fig.savefig("./density.pdf",format="pdf")
    

class null_distribution:

    def __init__(self):
        self.all_cps = np.array([])
        self.all_gc = np.array([])
        self.n = None
        self.sd = None
        self.bin_width = None
        self.mu_probDensByWidth = None
        self.ll_probDensByWidth = None
        self.initialized = False
        self.all_var = None
        self.all_mu = None
    
    def add_to_GC(self,gc_array,to_mask):
        gc_array = np.copy(gc_array) 
        
        for mask in to_mask:
            for i in xrange(np.shape(mask)[0]):
                gc_array[mask[i,0]:mask[i,1]] = -999.99
        
        where= gc_array!=-999.99
        gc_array = gc_array[where]
        self.all_gc = np.r_[self.all_gc,gc_array]

    def add(self,cp_array,to_mask):
        cp_array = np.copy(cp_array) 
        
        for mask in to_mask:
            for i in xrange(np.shape(mask)[0]):
                cp_array[mask[i,0]:mask[i,1]] = -999.99
        
        where= cp_array!=-999.99
        cp_array = cp_array[where]
        
        self.all_cps = np.r_[self.all_cps,cp_array]
        
        #sd = np.std(self.all_cps)
        #n = float(self.all_cps.shape[0])
        #bin_width = (3.5*sd)/(n**(1.0/3.0))
        #print "sd: %f, bin_width: %f, n_windows total %d"%(sd, bin_width,n)
        #bins = np.arange(-5*sd,5*sd,bin_width)
        #print "nbins:",bins.shape
        """
        1 - remove the gapped regions
        2 - determine the # of wnds??
        """

    def setup_mu_dist(self,width):
        """
        1d
        """
        last_chunk = self.csum_all[width:]
        first_chunk = self.csum_all[:-width]
        delta = (last_chunk-first_chunk)/width
        mu = np.mean(delta)
        sd = np.std(delta)
        n = delta.shape[0]
        bin_width = (3.5*sd)/(n**(1.0/3.0))
        bin_width /= 10 
        bin_edges = np.r_[np.array([-99999.9]), 
                          np.arange(-5*sd,5*sd,bin_width),   
                          np.array([99999.9])
                          ]

        density,edges = np.histogram(delta, bin_edges, density=True)
        
        density = density*np.diff(edges)
        cum_density= np.cumsum(density)

        self.mu_probDensByWidth[width] = tuple([cum_density, density, edges])

        #print "width %d, sd: %f n: %f binwidth: %f"%(width, sd, n, bin_width),np.diff(edges)[1]
    
    def _setup_mu_dist(self,width):
        """
        2d
        """
        
        last_chunk = self.csum_all[width:]
        first_chunk = self.csum_all[:-width]
        delta = (last_chunk-first_chunk)/width
        
        sd = np.std(delta)
        n = delta.shape[0]
        bin_width = (3.5*sd)/(n**(1.0/3.0))
        bin_width /= 10 
        bin_edges = np.r_[np.array([-99999.9]), 
                          np.arange(-5*sd,5*sd,bin_width),   
                          np.array([99999.9])
                          ]


        var_vect = (self.var_left+np.roll(self.var_right,-width))/2.0
        var_vect = var_vect[width:]
        #density,edges = np.histogram(delta, bin_edges, density=True)
        print bin_edges
        print self.variance_bins
        print delta.shape
        print var_vect.shape
        histogram, edges_mu, edges_var = np.histogram2d(delta, var_vect, bins=[bin_edges, self.variance_bins])
        print histogram.shape
        print np.sum(histogram,0) 
        density=histogram/np.sum(histogram,0) 
        print np.var(density,0)
        print density.sum(0)
        plot2d(density,edges_mu,edges_var)
        #raw_input()
        print width
        exit(1)

        density  = np.cumsum(density*np.diff(edges))
        self.mu_probDensByWidth[width] = tuple([density,edges])
        #print "width %d, sd: %f n: %f binwidth: %f"%(width, sd, n, bin_width),np.diff(edges)[1]
    
    def setup_ll_dist(self,width):
        
        last_chunk = self.csum_all_log_ps[width:]
        first_chunk = self.csum_all_log_ps[:-width]
        delta = (last_chunk-first_chunk)
        sd = np.std(delta)
        n = delta.shape[0]
        bin_width = (3.5*sd)/(n**(1.0/3.0))
        bin_width /= 10 
        bin_edges = np.r_[np.array([-99999.9]), 
                          np.arange(-5*sd,5*sd,bin_width),   
                          np.array([99999.9])
                          ]

        density,edges = np.histogram(delta, bin_edges, density=True)
        density  = np.cumsum(density*np.diff(edges))
        self.ll_probDensByWidth[width] = tuple([density,edges])
        #print "width %d, sd: %f n: %f binwidth: %f"%(width, sd, n, bin_width),np.diff(edges)[1]
    
    def get_indiv_window_dist(self):
        """
        get the distribution of single windows
        """
        mu = np.mean(self.all_cps)
        sd = np.std(self.all_cps)
        n = self.all_cps.shape[0]
        bin_width = (3.5*sd)/(n**(1.0/3.0))
        bin_width /= 10 
        bin_edges = np.r_[np.array([-99999.9]), 
                          np.arange(-5*sd,5*sd,bin_width),   
                          np.array([99999.9])
                          ]

        density,edges = np.histogram(self.all_cps, bin_edges, density=True)
        density  = np.cumsum(density*np.diff(edges))
        self.indiv_win_density = density
        self.indiv_win_edges = edges

        self.all_cp_ps  = density[np.clip(np.searchsorted(self.indiv_win_edges,self.all_cps)+1, 0, density.shape[0]-1)]
        self.all_cp_ps[self.all_cp_ps>0.5] = 1-self.all_cp_ps[self.all_cp_ps>0.5] 
        self.log_all_cp_ps = np.log(self.all_cp_ps)
        self.csum_all_log_ps = np.cumsum(self.log_all_cp_ps)
    
    
    def get_variance_dist(self):
        
        dec=10
        sorted_variance = np.sort(self.variance_vect)
        n = self.variance_vect.shape[0]
        n_dec = n/dec
        
        bins = []
        curr_dec = 0

        for i in xrange(dec):
            bins.append(sorted_variance[curr_dec])
            curr_dec+=n_dec

        bins.append(sorted_variance[-1])
        histogram, edges = np.histogram(self.variance_vect, bins)
        self.variance_bins = edges

        print 'variance decile edges:', edges

    def init_all_dists(self):
        
        self.get_indiv_window_dist()
        self.csum_all = np.cumsum(self.all_cps)
        self.mu_probDensByWidth = {}
        self.ll_probDensByWidth = {}

        self.half_width = 250
        self.variance_vect = get_windowed_variance(self.all_cps, self.half_width)
        self.var_left=np.roll(self.variance_vect,self.half_width+1)
        self.var_right=np.roll(self.variance_vect,-self.half_width-1)
        """
        var_left and var_right represent at position k, 
        the variance of the 2*half_width+1 windows to the right
        and the left 
        """
        self.get_variance_dist()
        self.null_var = np.var(self.all_cps)
        self.null_mu = np.mean(self.all_cps)
        self.initialized = True
    
    def get_ll_based_p_value(self,call):
        """
        get p_value based on the null distribution of 
        log-likelihood of consecutive windows with values
        """
         
        if not self.initialized:
            self.init_all_dists()
        
        
        width = call.wnd_end - call.wnd_start + 1
        values = call.values
        value = np.sum(np.log(self.indiv_win_density[ np.clip( np.searchsorted(self.indiv_win_edges, values)+1,
                                                 0,
                                                 self.indiv_win_density.shape[0]-1) 
                                                 ]))

        if not width in self.ll_probDensByWidth:
            t = time.time()
            self.setup_ll_dist(width)
        
        density, edges = self.ll_probDensByWidth[width]
        
        where = np.searchsorted(edges,value)+1
        where = min(where, density.shape[0] - 1)
        
        p = where < density.shape[0] and density[where] or 1.0
        
        #print value, width, p, call.start, call.end 
        #plot(density, edges, where) 
        #raw_input()
        
        if p>0.5: p=1-p
        if p<0: p=0
            
        return p
    

    def get_mu_based_p_value_simple(self, value, width):
        """
        get p_value based on the null distribution of means of 
        windows of same size
        """
         
        if not self.initialized:
            self.init_all_dists()
        
        #if width>10: return 0
        if not width in self.mu_probDensByWidth:
            self.setup_mu_dist(width)
        
        cum_density, density, edges = self.mu_probDensByWidth[width]
        
        where = np.searchsorted(edges,value)+1
        where = min(where, density.shape[0] - 1)
        
        p = where < density.shape[0] and cum_density[where] or 1.0
        
        #print value, width, p, call.start, call.end 
        #plot(density, edges, where) 
        #raw_input()
        
        if p>0.5: p=1-p
        if p<0: p=0
            
        return p

    def explore_mu_based_p_values(self):
        F = open('./distributions/calibrate.txt','w')
        print >>F, "width\tp\tmu"
        
        for width in xrange(1,100):
            if not width in self.mu_probDensByWidth:
                self.setup_mu_dist(width)
            
            cum_density, density, edges = self.mu_probDensByWidth[width]
            #plot_density(density, edges, "./distributions/%d.pdf"%width)
            for i in xrange(5,100,5):
                v=i/100.0
                p=self.get_mu_based_p_value_simple(v, width)
                print>>F, "%d\t%f\t%f"%(width, p, v)
    
    
    def explore_t_test_p_values(self):
        F = open('./distributions/calibrate.txt','w')
        print >>F, "width\tp\tmu"
        
        for width in xrange(1,100,10):
            
            for i in xrange(5,100,5):
                v=i/100.0
                p=self.get_t_test_p_value(v)
                print>>F, "%d\t%f\t%f"%(width, p, v)
    
    def get_rank_sum_p_value(self,call):
        z,p = stats.ranksums(call.values,self.all_cps)
        #return p
        return p
        
    def get_t_stats(self,call):
        
        if not self.initialized:
            self.init_all_dists()

        mu1 = np.mean(call.values)
        var1 = np.var(call.values)
        return {"null_var" : self.null_var, 
                "null_mu" : self.null_mu, 
                "null_n" : self.all_cps.shape[0],
                "mu1" : mu1, 
                "var1" : var1, 
                "n1" : call.values.shape[0] 
                }
    
    def get_all_var(self):
        if self.all_var == None:
            self.all_var = np.var(self.all_cps) 

        return self.all_var
    
    def get_all_mu(self):
        if self.all_mu == None:
            self.all_mu = np.mean(self.all_cps) 

        return self.all_mu


    def get_corrected_t_test_p_value(self,call):
        """
        http://en.wikipedia.org/wiki/Student's_t-test 
        students t test with unqequal var and sample size
        correcting t
        """
        #t,p = stats.ttest_ind(call.values, self.all_cps, equal_var=False)
        #n's are corrected by sqrt(n)
        #n1=call.values.shape[0]
        #n2=self.all_cps.shape[0]
        
        n1=np.sqrt(call.values.shape[0])
        n1=max(n1,2) #see the degrees of freedom, force it +ve
        n2=np.sqrt(self.all_cps.shape[0])
        
        var1=np.var(call.values)
        if var1 == 0.0: return 0.5
        #var2=np.var(self.all_cps)  #
        var2=self.get_all_var()

        mu1=np.mean(call.values)
        #mu2=np.mean(self.all_cps)  #
        mu2=self.get_all_mu()

        sX_X = np.sqrt( (var1/n1) + (var2/n2) )
        
        df_top = ((var1/n1) + (var2/n2))**2
        df_bottom = ((var1/n1)**2)/(n1-1) + ((var2/n2)**2)/(n2-1)

        df = df_top/df_bottom
        T = (mu1-mu2)/sX_X
        
        p = stats.t.cdf(T,df)
        p = min(p,1-p)
        #if p==0.0:  
        #    print p, n1, n2, var1, var2, mu1, mu2, sX_X, df, T
        return p 

    def get_t_test_p_value(self,call):
        """
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html#scipy.stats.ttest_ind
        """
        t,p = stats.ttest_ind(call.values, self.all_cps, equal_var=False)
        w = call.wnd_end - call.wnd_start +1
        return p
        

    def get_mu_based_p_value(self,call):
        """
        get p_value based on the null distribution of means of 
        windows of same size
        """
        width = call.wnd_end - call.wnd_start + 1
        value = call.value
        return min(self.get_mu_based_p_value_simple(value, width)*(width**0.5),0.5)
