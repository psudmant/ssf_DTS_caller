import numpy as np
from sys import stderr
from traverse_contours import get_cp_medians, get_cp_means


class edge:
    def __init__(self,pos,delta,left_epos,right_epos,left_median,right_median):
        self.pos = pos
        self.delta = delta
        self.left_epos = left_epos
        self.right_epos = right_epos
        self.left_median = left_median
        self.right_median = right_median
    
    def __str__(self):
        return "%d %f"%(self.pos,self.delta)
    
    def __repr__(self):
        return "%d %f l:%d r:%d"%(self.pos,self.delta,self.left_epos,self.right_epos)

def hierarch_merge_edges(copies,edges_passing_cutoff,max_merge_dif,use_mean=False):
    
    """
    note all points in the arrays here are just blocks, 
    they have no  start and end yet
    for n edges there are n-2 med_difs. Each med_dif 
    corresponds to an edge with the exception of the first 
    and last edge
    """
    np_copies=np.array(copies)
    np_edges_pass=np.array(edges_passing_cutoff)
    
    print "dtype"
    print np_copies
    print np_copies.dtype

    if use_mean:
        medians=get_cp_means(np_copies,np_edges_pass)
    else:
        medians=get_cp_medians(np_copies,np_edges_pass)
    print medians
    med_difs =np.abs(np.diff(medians))
    """ 
    edge at k corresponds to medians at k-1 and k
    median @ k corresponds to edges k and k+1
    median_dif @ l corresponds to edge at l+1 intersecting  medians  l and l+1 
    """
    arg_dif_sorted=np.argsort(med_difs)    

    ###First put the edges in order of into a list from smallest to largest delta
    edges_sorted_by_delta = []
    edge_pos_to_edge = {}
    #edge_pos_to_sorted_loc = {} #WON'T WORK< pos will change
    
    for i in xrange(arg_dif_sorted.shape[0]): ###ITERATING THROUGH THE DIFS
        delta_i_in_delta_list     = arg_dif_sorted[i]
        edge_i_in_edge_list     = arg_dif_sorted[i]+1
        delta=med_difs[delta_i_in_delta_list]

        edge_pos=edges_passing_cutoff[edge_i_in_edge_list]
        left_edge_pos =    edges_passing_cutoff[edge_i_in_edge_list-1]
        right_edge_pos =    edges_passing_cutoff[edge_i_in_edge_list+1]
        left_median=medians[edge_i_in_edge_list-1]
        right_median=medians[edge_i_in_edge_list]

        curr_edge=edge(edge_pos,delta,left_edge_pos,right_edge_pos,left_median,right_median)

        edges_sorted_by_delta.append(curr_edge)
        edge_pos_to_edge[edge_pos]=curr_edge
        #edge_pos_to_sorted_loc[edge_pos]=i

    #SPECIAL EDGES NOT TO BE SORTED
    edge_far_left=edge(edges_passing_cutoff[0],1000000,-1,edges_passing_cutoff[1],-1,medians[0])
    edge_pos_to_edge[edges_passing_cutoff[0]]=edge_far_left
    
    """
    added this to fix a bug in snowflake..
    not exactly sure why or how this works... 
    possiblly bad???, added Aug 7 2012.
    """

    edge_pos_to_edge[-1]=edge_far_left

    edge_far_right=edge(edges_passing_cutoff[-1],1000000,edges_passing_cutoff[-2],-1,medians[-1],-1) ##-1 was 0 for edges_passing_cutoff
    edge_pos_to_edge[edges_passing_cutoff[-1]]=edge_far_right
    ###ALSO, changking delta to 100000 so, it doesn't get merged    AND NOW, adding to the sorted edges
    edges_sorted_by_delta.append(edge_far_left)
    edges_sorted_by_delta.append(edge_far_right)

    edges_sorted_by_delta.reverse()
    """ 
    #NOW EDGE ARRAY IS SET UP IN ADDITION TO edge_pos_to_edge hash
    #ARRAY is reversed for fast end popping operations
    """ 
    min_dif = edges_sorted_by_delta[-1].delta
    stderr.write("iteration...")
    i=0
    while min_dif    < max_merge_dif:
        if i%100==0:
            stderr.write("%d.."%i)
            stderr.flush()
        i+=1
        """
        merge the segments intersected by min_dif
        pop out that edge
        """
        edge_to_merge = edges_sorted_by_delta.pop()
        left_epos=edge_to_merge.left_epos
        right_epos=edge_to_merge.right_epos
        
        left_edge=edge_pos_to_edge[left_epos]
        right_edge=edge_pos_to_edge[right_epos]

        #calc a new median
        if use_mean:
            new_median = np.mean(copies[left_epos:right_epos])
        else:
            new_median = np.median(copies[left_epos:right_epos])
        """
        update left and right edge median difs (deltas)
        update left and right edge left and right medians
        """

        left_edge.delta = abs(left_edge.left_median - new_median)
        left_edge.right_median = new_median

        right_edge.delta = abs(right_edge.right_median - new_median)
        right_edge.left_median = new_median
        
        #update left and right edge left and right e_poses
        right_edge.left_epos=left_epos
        left_edge.right_epos=right_epos

        #resort the array
        edges_sorted_by_delta = sorted(edges_sorted_by_delta,key=lambda x:-x.delta)
        min_dif = edges_sorted_by_delta[-1].delta
        
    print >>stderr,""
    #print edges_sorted_by_delta    

    final_edges=sorted(edges_sorted_by_delta,key=lambda x:x.pos)
    
    segments_s=[]
    segments_e=[]
    cps=[]
    #GC_info=defaultdict(list)
    
    for i in xrange(len(final_edges)-1):
        seg_start=final_edges[i].pos
        seg_end=final_edges[i+1].pos  #so- this is "i" because it will be accesssed int he ends, before it was accessed in the starts, so, it was the LAST start
        #seg_end=final_edges[i+1].pos #cojmmented May 13 2013
        segments_s.append(seg_start)
        segments_e.append(seg_end)
        cps.append(final_edges[i].right_median)
        #GC_info['max'].append(np.max(GC_content[seg_start:seg_end]))
        #GC_info['min'].append(np.max(GC_content[seg_start:seg_end]))
        #GC_info['med'].append(np.max(GC_content[seg_start:seg_end]))
        #GC_info['mean'].append(np.max(GC_content[seg_start:seg_end]))
    
    return segments_s,segments_e,cps
    

def hierarch_merge_edges_(copies,edges_passing_cutoff,max_merge_dif,use_mean=False):
    
    """
    DEPRECATED - only works with abutting edges
    for n edges there are n-2 med_difs. Each med_dif 
    corresponds to an edge with the exception of the first 
    and last edge
    """
    np_copies=np.array(copies)
    np_edges_pass=np.array(edges_passing_cutoff)
    if use_mean:
        medians=get_cp_means(np_copies,np_edges_pass)
    else:
        medians=get_cp_medians(np_copies,np_edges_pass)
    print medians
    med_difs =np.abs(np.diff(medians))
    """ 
    edge at k corresponds to medians at k-1 and k
    median @ k corresponds to edges k and k+1
    median_dif @ l corresponds to edge at l+1 intersecting  medians  l and l+1 
    """
    arg_dif_sorted=np.argsort(med_difs)    

    ###First put the edges in order of into a list from smallest to largest delta
    edges_sorted_by_delta = []
    edge_pos_to_edge = {}
    #edge_pos_to_sorted_loc = {} #WON'T WORK< pos will change
    
    for i in xrange(arg_dif_sorted.shape[0]): ###ITERATING THROUGH THE DIFS
        delta_i_in_delta_list     = arg_dif_sorted[i]
        edge_i_in_edge_list     = arg_dif_sorted[i]+1
        delta=med_difs[delta_i_in_delta_list]

        edge_pos=edges_passing_cutoff[edge_i_in_edge_list]
        left_edge_pos =    edges_passing_cutoff[edge_i_in_edge_list-1]
        right_edge_pos =    edges_passing_cutoff[edge_i_in_edge_list+1]
        left_median=medians[edge_i_in_edge_list-1]
        right_median=medians[edge_i_in_edge_list]

        curr_edge=edge(edge_pos,delta,left_edge_pos,right_edge_pos,left_median,right_median)

        edges_sorted_by_delta.append(curr_edge)
        edge_pos_to_edge[edge_pos]=curr_edge
        #edge_pos_to_sorted_loc[edge_pos]=i

    #SPECIAL EDGES NOT TO BE SORTED
    edge_far_left=edge(edges_passing_cutoff[0],1000000,-1,edges_passing_cutoff[1],-1,medians[0])
    edge_pos_to_edge[edges_passing_cutoff[0]]=edge_far_left
    
    """
    added this to fix a bug in snowflake..
    not exactly sure why or how this works... 
    possiblly bad???, added Aug 7 2012.
    """

    edge_pos_to_edge[-1]=edge_far_left

    edge_far_right=edge(edges_passing_cutoff[-1],1000000,edges_passing_cutoff[-2],-1,medians[-1],-1) ##-1 was 0 for edges_passing_cutoff
    edge_pos_to_edge[edges_passing_cutoff[-1]]=edge_far_right
    ###ALSO, changking delta to 100000 so, it doesn't get merged    AND NOW, adding to the sorted edges
    edges_sorted_by_delta.append(edge_far_left)
    edges_sorted_by_delta.append(edge_far_right)

    edges_sorted_by_delta.reverse()
    """ 
    #NOW EDGE ARRAY IS SET UP IN ADDITION TO edge_pos_to_edge hash
    #ARRAY is reversed for fast end popping operations
    """ 
    min_dif = edges_sorted_by_delta[-1].delta
    stderr.write("iteration...")
    i=0
    while min_dif    < max_merge_dif:
        if i%100==0:
            stderr.write("%d.."%i)
            stderr.flush()
        i+=1
        """
        merge the segments intersected by min_dif
        pop out that edge
        """
        edge_to_merge = edges_sorted_by_delta.pop()
        left_epos=edge_to_merge.left_epos
        right_epos=edge_to_merge.right_epos
        
        left_edge=edge_pos_to_edge[left_epos]
        right_edge=edge_pos_to_edge[right_epos]

        #calc a new median
        if use_mean:
            new_median = np.mean(copies[left_epos:right_epos])
        else:
            new_median = np.median(copies[left_epos:right_epos])
        """
        update left and right edge median difs (deltas)
        update left and right edge left and right medians
        """

        left_edge.delta = abs(left_edge.left_median - new_median)
        left_edge.right_median = new_median

        right_edge.delta = abs(right_edge.right_median - new_median)
        right_edge.left_median = new_median
        
        #update left and right edge left and right e_poses
        right_edge.left_epos=left_epos
        left_edge.right_epos=right_epos

        #resort the array
        edges_sorted_by_delta = sorted(edges_sorted_by_delta,key=lambda x:-x.delta)
        min_dif = edges_sorted_by_delta[-1].delta
        
    print >>stderr,""
    #print edges_sorted_by_delta    

    final_edges=sorted(edges_sorted_by_delta,key=lambda x:x.pos)
    
    segments_s=[]
    segments_e=[]
    cps=[]
    #GC_info=defaultdict(list)
    
    for i in xrange(len(final_edges)-1):
        seg_start=final_edges[i].pos
        seg_end=final_edges[i].pos  #so- this is "i" because it will be accesssed int he ends, before it was accessed in the starts, so, it was the LAST start
        #seg_end=final_edges[i+1].pos #cojmmented May 13 2013
        segments_s.append(seg_start)
        segments_e.append(seg_end)
        cps.append(final_edges[i].right_median)
        #GC_info['max'].append(np.max(GC_content[seg_start:seg_end]))
        #GC_info['min'].append(np.max(GC_content[seg_start:seg_end]))
        #GC_info['med'].append(np.max(GC_content[seg_start:seg_end]))
        #GC_info['mean'].append(np.max(GC_content[seg_start:seg_end]))
    
    return segments_s,segments_e,cps
    

def lr_merge_edges(copies,edges_passing_cutoff,max_merge_dif):
    """
    DEPRECATED
    """
    final_edges=[]
    edge_copies=[]

    final_edges.append(edges_passing_cutoff[0])
    final_edges.append(edges_passing_cutoff[1])
    edge_copies.append(np.median(copies[final_edges[0]:final_edges[-1]]))
    
    #copies has a large shape,
    #edges_passing is shorter, and refers to the windows in copies

    for i in xrange(2,len(edges_passing_cutoff)):
        curr_chunk_s = final_edges[-1]
        curr_chunk_e = edges_passing_cutoff[i]
        
        curr_chunk_cp = np.median(copies[curr_chunk_s:curr_chunk_e])
        if abs(edge_copies[-1] - curr_chunk_cp)>max_merge_dif:
            final_edges.append(curr_chunk_e)    
            edge_copies.append(curr_chunk_cp)
        else:
            final_edges[-1] = curr_chunk_e
            edge_copies[-1]=np.median(copies[final_edges[-2]:final_edges[-1]])
    
    segments_s=[]
    segments_e=[]
    cps=[]
    for i in xrange(len(final_edges)-1):
        segments_s.append(final_edges[i])
        segments_e.append(final_edges[i+1])
        cps.append(edge_copies[i])

    return segments_s,segments_e,cps

def get_GC_over_region(chr,start,end,fn_GC):
    
    n_expected=(end-start)/1000+2
    i=0
    gc_cont = np.zeros(n_expected,np.float)
    
    for line in tabix.Tabix(fn_GC).fetch(chr,start,end):

        chr,start,end,gc=line.split()
        gc_cont[i]=.1#$float(gc)    
        i+=1
    
    if i<n_expected:
        print >> stderr, "compressing GC ARRAY"
        new_gc_cont = np.zeros(n_expected-1,np.float)
        new_gc_cont=gc_cont[0:n_expected-1]
        gc_cont=new_gc_cont
        
    return np.median(gc_cont)

