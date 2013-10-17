import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mCols 
from matplotlib import cm as cm

import matplotlib.collections as collections
#from pytabix import tabix
import pysam
import numpy as np
from sys import stderr

class plot():
    def __init__(self,chr,caller,fn_gene_tabix,fn_dup_tabix,plot_lims,out_dir,callset):
        self.chr = chr
        self.caller = caller
        self.n_wavelet_scales=len(self.caller.scales)
        self.fn_gene_tabix = fn_gene_tabix
        self.fn_dup_tabix = fn_dup_tabix
        self.out_dir = out_dir
        self.plot_lims = plot_lims
        self.callset = callset
         
         
    def plot_all(self,chunk_len=1000000,bp_start=0,add_heatmap=False):
        total_chr_len = self.caller.ends[-1]
        bp_end = min(total_chr_len,bp_start+chunk_len)
        
        print "%s:0-%d"%(self.chr, total_chr_len)
        
        while bp_start<total_chr_len:
            
            print "%s:%d-%d"%(self.chr,bp_start, bp_end)
            #new_plot = line_plot(chr,bp_start,bp_end,self.caller,self.fn_gene_tabix,self.fn_dup_tabix)
            new_plot = line_plot(self.chr,bp_start,bp_end,self.caller,self.fn_gene_tabix,self.fn_dup_tabix,self.callset,add_heatmap=add_heatmap)
            new_plot.plot_features()
            new_plot.output("%s/%s_%d_%d"%(self.out_dir,self.chr,bp_start,bp_end),self.plot_lims)
            bp_start = bp_end
            bp_end = min(total_chr_len,bp_start+chunk_len)

class line_plot:
    
    def __init__(self,chr,start,end,caller,fn_gene_tabix,fn_dup_tabix,callset,add_heatmap=False):
        """
        .05 spacer
            80%
        .05 spacer
            20% 
        .05 spacer
        """
        ##### INIT PLOTTING PALLETTE ######
        total_annotations=2.0
        
        self.callset = callset

        self.fn_gene_tabix = fn_gene_tabix
        self.fn_dup_tabix = fn_dup_tabix
        self.height = 11
        self.width = 8.5
        self.add_heatmap=add_heatmap

        size = [self.width,self.height]
        plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
        
        self.fig = plt.figure()    
        self.fig.set_figwidth(size[0])
        self.fig.set_figheight(size[1])
        
        der_axes_frac=.4

        self.e_spacer=-.04
        self.caller = caller
        self.n_wavelet_scales=len(self.caller.scales)
        
        margin = 0.05 #of width
        width=1-margin*2
        
        main_plot_frac = .8
        annot_frac = 1-main_plot_frac
        total_spacers = 3*margin 
        
        t_non_spacer = 1-total_spacers
        
        indiv_annot_height = (annot_frac*t_non_spacer)/total_annotations
        cp_height  =  (main_plot_frac*t_non_spacer)*(1-der_axes_frac)
        der_height =  (main_plot_frac*t_non_spacer)*(der_axes_frac)

        curr_top  = 1-margin
        plot_rect=[margin,curr_top-cp_height,width,cp_height]
        curr_top-=cp_height
        der_plot_rect=[margin,curr_top-der_height,width,der_height]
        curr_top-=der_height+margin
            
        self.axes = self.fig.add_axes(plot_rect)
        self.der_axes = self.fig.add_axes(der_plot_rect)
        
        gene_rect=[margin,curr_top-indiv_annot_height,width,indiv_annot_height]
        self.gene_ax = self.fig.add_axes(gene_rect,frameon=False)
        curr_top-=indiv_annot_height
        dup_rect=[margin,curr_top-indiv_annot_height,width,indiv_annot_height]
        self.dup_ax = self.fig.add_axes(dup_rect,frameon=False)
        
        if self.add_heatmap:
            self.heatmap_fig = plt.figure()    
            self.heatmap_fig.set_figwidth(size[0])
            self.heatmap_fig.set_figheight(size[1])
            
            sep = .1
            height=(1-(3*sep))/2.0
            
            plot_rect_1=[margin,1-sep-height,width,height]
            plot_rect_2=[margin,1-(2*sep)-(2*height),width,height]
            
            self.heat_ax = self.heatmap_fig.add_axes(plot_rect_1,frameon=False)
            self.heat_ax_dir = self.heatmap_fig.add_axes(plot_rect_2,frameon=False)

        ####### INIT REGION COORDINATES ######
        self.chr = chr
        self.start = start
        self.end = end
    
    def plot_features(self):
        self.plot_dups(self.fn_dup_tabix)
        self.plot_genes(self.fn_gene_tabix)
        
            
        wnd_start = np.searchsorted(self.caller.starts,self.start)
        wnd_end = np.searchsorted(self.caller.starts,self.end)
        if self.add_heatmap:
            self.plot_magnitude_heatmap(wnd_start,wnd_end)
        
        self.plot_wavelet_convolution(wnd_start,wnd_end)

        curr_scale_k=0
        
        #for scale in caller.scales:
        #    (edges,pos_edges,neg_edges)=caller.transitions_by_wnd_scale[wnd_size][scale]
        #    new_plot.add_wavelet_edge(k,curr_scale_k,wnd_start,wnd_end,starts,scale,edges,pos_edges,neg_edges)
        #    curr_scale_k+=1    

    def plot_dups(self,fn_dups):
            
        #dups = regions_from_bed(fn_dups)
        #locs,vals = dups.get_locations_over_interval(self.chr,self.start,self.end)
        #line_segs = collections.LineCollection(locs)

        max_loc = -1
        y=0
        print self.chr,self.start,self.end
        #for dup_line in tabix.Tabix(fn_dups).fetch(self.chr,self.start,self.end):
        #for i in xrange(locs.shape[0]):
        tbx_dups = pysam.Tabixfile(fn_dups) 
        for dup_line in tbx_dups.fetch(self.chr,self.start,self.end):
            sline=dup_line.split()
            pc_id=float(sline[26])
            start,end=int(sline[2]),int(sline[3])
            y=start<max_loc and y-.5 or 0
            max_loc = max(max_loc,end)
            if pc_id>0.93:
                self.dup_ax.plot([start,end],[y,y],'orange',linewidth=1)    
            else:
                self.dup_ax.plot([start,end],[y,y],'grey',linewidth=1)    
    
    def plot_genes(self,fn_genes):
            
        #genes = regions_from_bed(fn_genes,names=True)
        #locs,vals = genes.get_locations_over_interval(self.chr,self.start,self.end)
        
        max_loc = -1
        y=0
        last_start,last_end=0,0
        #for i in xrange(locs.shape[0]):

        count=0
        max_count=10
        #for gene_line in tabix.Tabix(fn_genes).fetch(self.chr,self.start,self.end):
        tbx_genes = pysam.Tabixfile(fn_genes) 
        for gene_line in tbx_genes.fetch(self.chr,self.start,self.end):
            sline=gene_line.split()
            start,end,name=int(sline[4]),int(sline[5]),sline[12]
            if start==last_start and end==last_end: continue
            #y=start<max_loc and y-1.5 or 0
            if count==max_count:
                y=0
                count=0
            else:    
                y-=1
                count+=1
            
            f_size=8
            max_loc = max(max_loc,end+f_size*1.1*len(name))
            self.gene_ax.plot([start,end],[y,y],'g',linewidth=4,alpha=.6)    
            self.gene_ax.annotate(name,(end+100,y),fontsize=f_size,horizontalalignment='left')

            last_start,last_end=start,end

    def output(self,outfile,plot_lims):
    
        start = self.start
        end = self.end
        
        if plot_lims!=None:
            ylim_min,ylim_max = [float(x) for x in plot_lims.split(":")]
        else:    
            ylim_max=30
            ylim_min=0


        for k in range(int(ylim_max)):
            self.axes.plot(np.array([start,end]),np.array([k,k]),linewidth=.1,alpha=1,color='k')
            self.axes.plot(np.array([start,end]),np.array([-k,-k]),linewidth=.1,alpha=1,color='k')
        
        self.axes.set_xlim(start,end)
        self.axes.get_xaxis().set_ticks([])
        self.der_axes.set_xlim(start,end)
        self.axes.set_ylim(ylim_min,ylim_max)
        self.der_axes.set_ylim(self.e_spacer*self.n_wavelet_scales,-1*self.e_spacer)
        #self.der_axes.set_yscale('log')
        
        self.dup_ax.set_xlim(start,end)
        self.gene_ax.set_xlim(start,end)
        
        dup_ax_y_lim = self.dup_ax.get_ylim()
        self.dup_ax.set_ylim(dup_ax_y_lim[0]-.3,dup_ax_y_lim[1]+.3)
        
        self.fig.savefig("%s.pdf"%(outfile.rstrip(".pdf")),format='pdf',height=self.height,width=self.width)
        plt.close(1)     
        
        if self.add_heatmap:
            #self.heat_ax.set_xlim(start,end)
            self.heatmap_fig.savefig("%s_HM.png"%(outfile.rstrip(".pdf")),format='png',height=self.height,width=self.width)
            
    def add_wavelet_edge(self,k,curr_scale_k,wnd_start,wnd_end,x,scale,edges,pos_edges,neg_edges):
        #self.axes[k].plot(x[wnd_start:wnd_end],(edges[wnd_start:wnd_end]*self.e_spacer)+curr_scale_k*self.e_spacer,alpha=.3,color='k')
        ###POS EDGES ARE RED
        xs=x[wnd_start:wnd_end]
        p_edges=pos_edges[wnd_start:wnd_end]
        n_edges=neg_edges[wnd_start:wnd_end]

        p_xs = xs[np.where(p_edges!=0)]
        n_xs = xs[np.where(n_edges!=0)]
        ps = p_edges[np.where(p_edges!=0)]
        ns = n_edges[np.where(n_edges!=0)]

        self.der_axes[k].plot(x[wnd_start:wnd_end],(pos_edges[wnd_start:wnd_end]*self.e_spacer)+curr_scale_k*self.e_spacer,alpha=.5,color='r',linewidth=.05)
        self.der_axes[k].plot(p_xs,(ps*self.e_spacer)+curr_scale_k*self.e_spacer,color='r',marker='.',alpha=.5,linewidth=.05,ms=.5,ls='None',mfc='r',mew=.1)
        ####NEG_EDGES ARE GREEN
        self.der_axes[k].plot(x[wnd_start:wnd_end],(neg_edges[wnd_start:wnd_end]*self.e_spacer)+curr_scale_k*self.e_spacer,alpha=.5,color='g',linewidth=.05)
        self.der_axes[k].plot(n_xs,(ns*self.e_spacer)+curr_scale_k*self.e_spacer,color='g',marker='.',alpha=.5,linewidth=.05,ms=.5,ls='None',mfc='g',mew=.1)
    #def plot_wavelet_convolution(self,k,wnd_start,wnd_end,x,g1,g2,cp,edges,smoothed_cp):
        
    #def plot_wavelet_convolution(self,k,wnd_start,wnd_end,x,cp,smoothed_cp,contour_intersects,cutoff_scale,segs_s,segs_e,seg_cps):
    
    def plot_magnitude_heatmap(self,wnd_start,wnd_end):

        starts = self.caller.starts[wnd_start:wnd_end]
        ends = self.caller.ends[wnd_start:wnd_end]

        ders = self.caller.der1[:,wnd_start:wnd_end]
        magnitudes = ders.ravel() 
        print "RANGE:",np.min(ders), np.max(ders) 
        dir_ders = ders.copy()  
        dir_ders[dir_ders>0] = 1 
        dir_ders[dir_ders<0] = -1 

        n_scales = self.caller.der1.shape[0]
        n_xs = starts.shape[0]  
        
        xs = np.tile((starts+ends)/2.0,n_scales)
        ys = np.repeat(np.arange(n_scales),n_xs)
        
        print xs, ys, magnitudes
        print "pcoloring" 
        self.heat_ax.pcolor(ders,cmap=cm.RdBu,vmax=.1, vmin=-.1)
        self.heat_ax_dir.pcolor(dir_ders,cmap=cm.RdBu,vmax=2,vmin=-2)
        print "done" 
        
        #self.hist2d.pcolormesh(xs,ys,magnitudes, shading='gouraud')
        #gridsize=n_xs
        #self.heat_ax.hexbin(xs,ys,C=magnitudes,gridsize=gridsize, cmap=CM.jet, bins=None)
        #self.heat_ax.hexbin(xs,ys,C=magnitudes, cmap=CM.jet, bins=None)
        
            

    def plot_wavelet_convolution(self,wnd_start,wnd_end):
        starts = self.caller.starts
        ends = self.caller.ends
        
        all_x_coords = (starts+ends)/2.0
        x_coords = (starts[wnd_start:wnd_end] + ends[wnd_start:wnd_end])/2.0
        self.axes.plot(x_coords,self.caller.cp_data[wnd_start:wnd_end],alpha=.8,color='b',linewidth=.3)
        self.axes.plot(x_coords,self.caller.cp_data[wnd_start:wnd_end],alpha=.8,color='b',linewidth=.3,marker='o', markersize=1)
        
        if x_coords.shape[0] ==0:
            print >>stderr, "BAILING OUT - in a gap"
            return

        self.der_axes.plot([x_coords[0],x_coords[-1]],[self.caller.cutoff_scale*self.e_spacer,self.caller.cutoff_scale*self.e_spacer],alpha=.5,color='b')
        
        for scale,intersects in self.caller.contour_intersects.iteritems():
            intersects =np.array(intersects)
            intersect_in_range=intersects[np.where((intersects>=wnd_start)&(intersects<wnd_end))[0]]
            l=intersect_in_range.shape[0]
            if l==0: continue
            xs=np.reshape(np.c_[all_x_coords[intersect_in_range],all_x_coords[intersect_in_range],all_x_coords[intersect_in_range]],(1,-1))[0]
            ys=np.reshape(np.c_[np.zeros(l),np.ones(l)*self.e_spacer*scale,np.zeros(l)],(1,-1))[0]
            #print xs,ys
            self.der_axes.plot(xs,ys,alpha=.8,color='k',linewidth=.1)

        segs_s, segs_e, seg_cps = self.caller.segment_edges 
        segs_s,segs_e = np.array(segs_s),np.array(segs_e)
        locs=np.where( ((segs_s<wnd_start)&(segs_e>wnd_start))|
                                        ((segs_s<wnd_end)&(segs_e>wnd_end))|
                                        ((segs_s>=wnd_start)&(segs_e<=wnd_end)) )[0]
        print locs
        
        if locs.shape[0]==0: 
            print >>stderr, "BAILING OUT - locs.shape is empty"
            return
        
        mn=np.amin(locs)
        mx=np.amax(locs)+1
        
        mx=min(mx,(len(segs_s)-1)) #in the last window we fail unless I put this
        
        starts_in_range=segs_s[mn:mx]
        ends_in_range=segs_e[mn:mx]
        cps_in_range=seg_cps[mn:mx]
        #starts_in_range=segs_s[locs]
        #ends_in_range=segs_e[locs]
        #cps_in_range=seg_cps[locs]
        l=starts_in_range.shape[0]
        
        print "starts",starts_in_range
        print "ends",ends_in_range
        
        #xs=[self.caller.starts[starts_in_range],self.caller.ends[ends_in_range]]
        #xs=[self.caller.starts[starts_in_range],self.caller.starts[ends_in_range]]
        xs=[all_x_coords[starts_in_range],all_x_coords[ends_in_range]]
        ys=[cps_in_range,cps_in_range]
        self.axes.plot(xs,ys,alpha=.8,color='r',linewidth=.5)
         
        for i,start in enumerate(starts[starts_in_range]): 
            self.axes.annotate(start,xy=(start,.5),xytext=(start,1),rotation=45,arrowprops=dict(facecolor='black',width=.002,headwidth=.005,frac=.1),fontsize=3)
            #self.gene_ax.annotate(name,(end+100,y),fontsize=f_size,horizontalalignment='left')
        
        calls  = self.callset.get_calls_in_range(wnd_start,wnd_end)
        for call in calls:
            x=(call.start+call.end)/2.0
            if call.fdr_significant:
                self.axes.annotate("%e.5 %d"%(call.p_value,call.significance_level), xy=(x,-.5), xytext=(x,-1), rotation=45, ha='right', 
                                   arrowprops=dict(color='red', facecolor='red',width=.002,headwidth=1,frac=.5),
                                   fontsize=5,
                                   color='red')
            else:
                self.axes.annotate("%e.5 %d"%(call.p_value,call.significance_level), xy=(x,-.5), xytext=(x,-1), rotation=45, ha='right', 
                                   arrowprops=dict(color='blue',facecolor='blue',width=.002,headwidth=1,frac=.5),
                                   fontsize=3,
                                   color='blue')

    def plot_CN_summary(self,wnd_DTS_by_genome,genome_to_plot_group,plot_groups_to_indivs,wnd_start,wnd_end,starts,color_hash):

        for plot_group,indivs in plot_groups_to_indivs.iteritems():
            n_indivs=len(indivs)
            cp_stack = np.zeros([n_indivs,wnd_end-wnd_start])
            c=color_hash[plot_group.upper()]
            stderr.write(plot_group)
            for i, indiv in enumerate(indivs):
                cp_stack[i,:]=wnd_DTS_by_genome[indiv]['copy'][self.chr][wnd_start:wnd_end]
                stderr.write(".")
                stderr.flush()
            print >>stderr,""
            mu=np.mean(cp_stack,0)
            sd=np.std(cp_stack,0)
            mn=np.min(cp_stack,0)
            mx=np.max(cp_stack,0)
            self.cp_ax.plot(starts[wnd_start:wnd_end],mu,color=c)
            self.cp_ax.fill_between(starts[wnd_start:wnd_end],mu-sd,mu+sd,alpha=.2,color=c)
            self.cp_ax.plot(starts[wnd_start:wnd_end],mn,linestyle='dot',alpha=1,color=c,size=.1)
        #FOR indiv in each group, make a stack
        #plot the mean, sd shade, min max dotted
        return 0    

    def plot_legend(self,color_hash):
        #self.cp_ax.legend()            
    
        pops = []
        lines= []

        for plot_group, color in color_hash.iteritems(): 
            pops.append(plot_group)
            lines.append(mpl.lines.Line2D([0,1],[0,0],color=color,linewidth=6))
        self.cp_ax.legend(lines,pops,ncol=4,loc=2,mode='expand',prop=mpl.font_manager.FontProperties(size=6))

