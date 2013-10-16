
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mCols 
import matplotlib.collections as collections
import numpy as np

class line_plot:
    
    def __init__(self,chr,start,end,n_plots,n_scales):

        ##### INIT PLOTTING PALLETTE ######
        total_annotations=2
        self.height = 5*(n_plots+total_annotations)
        self.width = 8.5
        size = [self.width,self.height]
        self.n_plots=n_plots    
        plt.rc('grid',color='0.75',linestyle='l',linewidth='0.1')
        self.fig = plt.figure()    
         self.fig.set_figwidth(size[0])
         self.fig.set_figheight(size[1])
        axescolor  = '#f6f6f6'
        

        der_axes_frac=.4

        self.e_spacer=-.04
        self.n_wavelet_scales=n_scales
        h_margin = 0.1 #of width
        
        #w_v_margin = 0.1 #of one window
        spacer_frac=.1 #of one window
        
        width=1-h_margin*2

        #extra_spaces=1
        t_height=1.0/(n_plots+(.5*total_annotations))
    
        height=(1-spacer_frac)*t_height
        spacer=spacer_frac*t_height

        #height=(1-(v_margin*2-spacer*(total_annotations+extra_spaces))-spacer*n_plots)/n_plots
        
        self.axes=[]
        self.der_axes=[]
        for i in xrange(n_plots):
            #plot_rect=[h_margin,1-((i+1)*(height+spacer)),width,height] pre having der axes
            plot_rect=[h_margin,1-((i+1)*(height+spacer))+(height*der_axes_frac),width,height*(1-der_axes_frac)]
            der_plot_rect=[h_margin,1-((i+1)*(height+spacer)),width,height*(der_axes_frac)]
            
            self.axes.append(self.fig.add_axes(plot_rect))
            self.der_axes.append(self.fig.add_axes(der_plot_rect))
        
        last_axis_h=1-(n_plots*(height+spacer))

        gene_rect=[h_margin,last_axis_h-(height/3+spacer),width,height/3.0]
         self.gene_ax = self.fig.add_axes(gene_rect,frameon=False)
        dup_rect=[h_margin,last_axis_h-2.0*(height/3+spacer),width,height/3.0]
         self.dup_ax = self.fig.add_axes(dup_rect,frameon=False)

        #plot_hfrac=.6
        #plot_height=plot_hfrac*height
        #plot_bottom=1-v_margin-plot_height

        #annot_hfrac=(1-plot_hfrac)/total_annotations
        #annot_height=annot_hfrac*height
        #gene_bottom=1-v_margin-spacer-spacer*extra_spaces-plot_height-annot_height
        #dup_bottom=1-v_margin-2*spacer-spacer*extra_spaces-plot_height-2*annot_height

        #plot_rect=[h_margin,plot_bottom,width,plot_height]
        #gene_rect=[h_margin,gene_bottom,width,annot_height]
        #dup_rect=[h_margin,dup_bottom,width,annot_height]

        #self.cp_ax = self.fig.add_axes(plot_rect)
         #self.dup_ax = self.fig.add_axes(dup_rect,frameon=False)
         
        ####### INIT REGION COORDINATES ######
        self.chr = chr
        self.start = start
        self.end = end

    def plot_dups(self,fn_dups):
            
        #dups = regions_from_bed(fn_dups)
        #locs,vals = dups.get_locations_over_interval(self.chr,self.start,self.end)
        #line_segs = collections.LineCollection(locs)

        max_loc = -1
        y=0
        print self.chr,self.start,self.end
        for dup_line in tabix.Tabix(fn_dups).fetch(self.chr,self.start,self.end):
        #for i in xrange(locs.shape[0]):
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
        for gene_line in tabix.Tabix(fn_genes).fetch(self.chr,self.start,self.end):
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


        for i in xrange(self.n_plots):
            for k in range(int(ylim_max)):
                self.axes[i].plot(np.array([start,end]),np.array([k,k]),linewidth=.1,alpha=1,color='k')
            for k in range(0,int(ylim_min),-1):
                self.axes[i].plot(np.array([start,end]),np.array([k,k]),linewidth=.1,alpha=1,color='k')
            self.axes[i].set_xlim(start,end)
            self.axes[i].get_xaxis().set_ticks([])
            self.der_axes[i].set_xlim(start,end)
            #self.axes[i].set_ylim(self.e_spacer*self.n_wavelet_scales-2,ylim_max)
            self.axes[i].set_ylim(ylim_min,ylim_max)
            self.der_axes[i].set_ylim(self.e_spacer*self.n_wavelet_scales,-1*self.e_spacer)
        
        self.dup_ax.set_xlim(start,end)
        self.gene_ax.set_xlim(start,end)
        
        dup_ax_y_lim = self.dup_ax.get_ylim()
        self.dup_ax.set_ylim(dup_ax_y_lim[0]-.3,dup_ax_y_lim[1]+.3)
        
        self.fig.savefig("%s.pdf"%(outfile.rstrip(".pdf")),format='pdf',height=self.height,width=self.width)
        plt.close(1)     
        #plot_hfrac=.6
        #plot_height=plot_hfrac*height
        #plot_bottom=1-v_margin-plot_height
        #self.cp_ax.plot([start,end],[2.0,2.0],'r--')
        #self.cp_ax.set_xlim(start,end)

        #self.cp_ax.set_ylim(0,min(self.cp_ax.get_ylim()[1],50))
        #
        #self.gene_ax.set_xlim(start,end)
        #self.dup_ax.set_xlim(start,end)

        #######adjust the y limits of the gene
        #######and dup tracks to make it more visible
        #gene_ax_y_lim = self.gene_ax.get_ylim()
        #dup_ax_y_lim = self.dup_ax.get_ylim()
        ##other_ax_y_lim = self.other_ax.get_ylim()
        #self.gene_ax.set_ylim(gene_ax_y_lim[0]-2,gene_ax_y_lim[1]+1)
        #self.dup_ax.set_ylim(dup_ax_y_lim[0]-2,dup_ax_y_lim[1]+1)

        ##self.other_ax.set_ylim(other_ax_y_lim[0]-1,other_ax_y_lim[1]+10)
    
        ######remove all ticks
        #self.dup_ax.get_xaxis().set_ticks([])
        #self.dup_ax.get_yaxis().set_ticks([])
        #self.gene_ax.get_xaxis().set_ticks([])
        #self.gene_ax.get_yaxis().set_ticks([])
        #
        ##curr_x_ax.set_ticks([])    
        ##self.other_axis.get_xaxis().set_ticks([])
            
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
        
    def plot_wavelet_convolution(self,k,wnd_start,wnd_end,x,cp,smoothed_cp,contour_intersects,cutoff_scale,segs_s,segs_e,seg_cps):

        print "plotting", k    
        #print x[wnd_start:wnd_end],cp[wnd_start:wnd_end]
        self.axes[k].plot(x[wnd_start:wnd_end],smoothed_cp[wnd_start:wnd_end],alpha=.5,color='m')
        self.axes[k].plot(x[wnd_start:wnd_end],cp[wnd_start:wnd_end],alpha=.8,color='b',linewidth=.3)
        self.axes[k].plot(x[wnd_start:wnd_end],cp[wnd_start:wnd_end],alpha=.8,color='b',marker='.',ms=.4)
        #self.axes[k].plot(x[wnd_start:wnd_end],g1[wnd_start:wnd_end],alpha=.8,color='r')
        #self.axes[k].plot(x[wnd_start:wnd_end],g2[wnd_start:wnd_end],alpha=.8,color='g')
        #self.axes[k].stem(x[edges],np.ones(edges.shape[0])*-10)
        self.der_axes[k].plot([x[wnd_start],x[wnd_end]],[cutoff_scale*self.e_spacer,cutoff_scale*self.e_spacer],alpha=.5,color='b')
        for scale,intersects in contour_intersects.iteritems():
            intersects =np.array(intersects)
            intersect_in_range=intersects[np.where((intersects>=wnd_start)&(intersects<wnd_end))[0]]
            l=intersect_in_range.shape[0]
            if l==0: continue
            xs=np.reshape(np.c_[x[intersect_in_range],x[intersect_in_range],x[intersect_in_range]],(1,-1))[0]
            ys=np.reshape(np.c_[np.zeros(l),np.ones(l)*self.e_spacer*scale,np.zeros(l)],(1,-1))[0]
            #print xs,ys
            self.der_axes[k].plot(xs,ys,alpha=.8,color='k',linewidth=.1)

        #segs_s,segs_e,seg_cps
        #final_edges=np.array(final_edges)
        #final_copies=np.array(final_copies)
    
        segs_s,segs_e = np.array(segs_s),np.array(segs_e)
        locs=np.where( ((segs_s<wnd_start)&(segs_e>wnd_start))|
                                        ((segs_s<wnd_end)&(segs_e>wnd_end))|
                                        ((segs_s>=wnd_start)&(segs_e<=wnd_end)) )[0]
        print locs
        #starts=np.where((segs_s<wnd_start)&(segs_e>wnd_end))[0]
        #ends=np.where((segs_e>=wnd_start)&(segs_e<wnd_end))[0]
        #if starts.shape==0 and ends.shape==0: return    
        if locs.shape[0]==0: return
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
        print len(x)
        for i in np.arange(l):
            xs=[x[starts_in_range[i]],x[ends_in_range[i]]]
            ys=[cps_in_range[i],cps_in_range[i]]
            ys2=[0.1,0.1]
            print xs,ys
            self.axes[k].plot(xs,ys,alpha=.8,color='r',linewidth=.5)
            self.axes[k].plot(xs,ys2,alpha=.8,color='g',linewidth=1)

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
            self.cp_ax.plot(starts[wnd_start:wnd_end],mn,linestyle='dot',alpha=1,color=c)
            self.cp_ax.plot(starts[wnd_start:wnd_end],mx,linestyle='dot',alpha=1,color=c)
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
