# distutils: language = c++
# distutils: sources = ./heap/heap.cc  ./edge.cc

import numpy as np
from sys import stderr
from traverse_contours import get_cp_medians, get_cp_means
import time

cdef extern from "heap.h":
    cdef cppclass heap_element[T]:

        heap_element(T)
        heap_element()
        
        T data
        unsigned int node_idx
         
        bint operator<(heap_element)

    cdef cppclass heap[T]:
        
        heap()
        #heap *new_heap "new heap" ()
        #void del_heap "delete" (heap *myheap)
        void init_heap(int)
        
        int insert(heap_element *)
        int remove(int)
        int remove(heap_element *)
        
        void print_v(int, int, int)
        void print_h(int)
        
        T peek_element(int)
        heap_element * get_parent(heap_element *)
        heap_element * get_left_child(heap_element *)
        heap_element * get_right_child(heap_element *)
        heap_element * get_element(int)
        
        int get_parent(int)
        int get_left_child(int)
        int get_right_child(int)
        int heap_integrity()
        
cdef extern from "edge.h":
    
    cdef cppclass edge(heap):

        edge()
        edge(long int, long int, long int, float, float, float)
        void edge_init(long int, long int, long int, float, float, float)
        void recalculate_delta()
        long int pos
        long int left_epos
        long int right_epos
        float delta
        float left_median
        float right_median
        edge *left_edge
        edge *right_edge
        
"""
ASSUME EDGES ARE >> inner THEN outer<< EDGES
thus, the if edges i and i +1 span
windows starting at k, k+1, k+2...k+l, 
then, copy for seg i->i+1
should be mean(k:k+l), ie, not including the 
k +lth wnd

so,  if each point is a window
        ._._._.
_._._._/       \._._._.__
        *       *
'.' indicates edge

this is how the get_cp_means works now
"""
    
def c_hierarch_merge_edges(cp_data,
                           edges_passing_cutoff,
                           max_merge,
                           use_mean,
                           max_pos,
                           starts,
                           ends):
    
    np_edges_pass = np.array(edges_passing_cutoff)
    np_copies=np.array(cp_data)

    """
    medians holds the mean/median value of block
    between edges, inclusive, exclusive, ie, 
    inner then outer coords 
    """

    medians = None
    
    if use_mean:
        medians=get_cp_means(np_copies,np_edges_pass)
    else:
        medians=get_cp_medians(np_copies,np_edges_pass)
   
    med_difs = np.abs(np.diff(medians))

    """ 
    edge at k corresponds to medians at k-1 and k
    median @ k corresponds to edges k and k+1
    median_dif @ l corresponds to edge at l+1 intersecting  medians  l-1 and l+1 
    """
    
    cdef heap[float] edge_heap
    edge_heap.init_heap(len(edges_passing_cutoff)+5)
    
    cdef heap_element[float]* h_elem
    cdef edge *last_edge
    cdef edge *curr_edge
    cdef edge *begin_edge
    cdef edge *end_edge
    cdef int i
    
    t=time.time()
    
    # create a fake edge on the far left
    curr_edge = new edge(0,-1,-1,9e9,9e9,medians[0])
    curr_edge.left_edge = NULL
    edge_heap.insert(<heap_element[float]*>curr_edge)
    last_edge = curr_edge
    begin_edge = curr_edge

    for i in xrange(len(edges_passing_cutoff)-1):
        edge_pos = edges_passing_cutoff[i]
        delta = i==0 and 9e9 or med_difs[i-1]
        curr_edge = new edge(edge_pos, -1, -1, delta, medians[i-1], medians[i])
        
        curr_edge.left_edge = last_edge
        last_edge.right_edge = curr_edge
        edge_heap.insert(<heap_element[float]*>curr_edge)
        last_edge = curr_edge
    
    # create a fake edge on the far right
    curr_edge = new edge(max_pos,-1,-1,9e9,medians[-1],9e9)
    curr_edge.left_edge = last_edge
    curr_edge.right_edge = NULL
    edge_heap.insert(<heap_element[float]*>curr_edge)
    last_edge.right_edge = curr_edge
    end_edge = curr_edge

    print "loading time: %fs"%(time.time()-t) 

    #now iterate over 'all nodes'
    cdef float new_med
    cdef int left_epos, right_epos 
    t = time.time() 
    for i in range( len( edges_passing_cutoff ) ):
        if i%100==0:
            stderr.write("%d.."%i)
            stderr.flush()
        
        if edge_heap.peek_element(0)>=max_merge:
            break
        
        curr_edge = <edge*>edge_heap.get_element(0)
        left_epos = curr_edge.left_edge.pos
        right_epos = curr_edge.right_edge.pos
        
        #get new median between left edge and right edge
        if use_mean:
            new_med = np.mean(np_copies[left_epos:right_epos])
        else:
            new_med = np.median(np_copies[left_epos:right_epos])
        
        ret = edge_heap.remove(<heap_element[float]*>curr_edge)
        ret = edge_heap.remove(<heap_element[float]*>curr_edge.right_edge)
        ret = edge_heap.remove(<heap_element[float]*>curr_edge.left_edge)
        
        curr_edge.left_edge.right_median = new_med
        curr_edge.left_edge.recalculate_delta()
        curr_edge.right_edge.left_median = new_med
        curr_edge.right_edge.recalculate_delta()

        curr_edge.right_edge.left_edge = curr_edge.left_edge
        curr_edge.left_edge.right_edge = curr_edge.right_edge
        
        edge_heap.insert(<heap_element[float]*>curr_edge.right_edge)
        edge_heap.insert(<heap_element[float]*>curr_edge.left_edge)
        
    stderr.write("\n")
    print "total merging time %fs"%(time.time()-t)
    
    #return all edges (start,end pairs) and their medians
    
    segment_starts = [] 
    segment_ends = [] 
    segment_cps = []
    
    curr_edge = begin_edge
    while curr_edge != end_edge:
        
        segment_starts.append(curr_edge.pos)
        segment_ends.append(curr_edge.right_edge.pos)
        segment_cps.append(curr_edge.right_median)
        curr_edge = curr_edge.right_edge

    return segment_starts, segment_ends, segment_cps


