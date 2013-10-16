# distutils: language = c

import numpy as np
cimport numpy as np
#import sys
#cimport cython 

#@cython.boundscheck(False)
#@cython.wraparound(False)

cdef int trace_contour(int scale_start,int x_start,np.ndarray[np.int8_t,ndim=2] d2_ints):
   
   cdef int curr_scale=scale_start      
   cdef int curr_sign=d2_ints[scale_start,x_start]
   
   cdef int search_i=0
   
   cdef int curr_x=x_start

   while curr_scale!=0:
      
      d2_ints[curr_scale,curr_x]=0 #zero it out. We've seen it
      ###LOOK below,left,right could optimize this with an array that I sort
      search_i=0
      while 1:
         #if search_i>curr_x or search_i+curr_x>d2_ints.shape[1]:
         #   break
         if curr_x-search_i==0: return 0
         if curr_x+search_i==d2_ints.shape[1]-1: return d2_ints.shape[1]-1

         if d2_ints[curr_scale-1,curr_x-search_i]==curr_sign:
            curr_x=curr_x-search_i   
            break
         elif d2_ints[curr_scale-1,curr_x+search_i]==curr_sign:
            curr_x=curr_x+search_i   
            break
         else:
            search_i+=1
      
      curr_scale-=1
         
   return curr_x
      

def get_contours(np.ndarray[np.int8_t,ndim=2] d2_ints):
    
    
    """
    returns:
        rets
            a dictionary keyed off of the scale
            containing a list for each scale of the x
            intercepts traced down to
        x_intersect_to_scale
            each intersect to it's scale
            could use this to prioritize the intercepts
    """
    #cdef Py_ssize_t i THIS IS WHAT I USED TO USE FOR MY ITERATORS 
    #apparently this is only for "purists"

    rets = {}
    x_intersect_to_scale = {}
    cdef int curr_scale
    cdef int curr_x_in_scale
    for curr_scale in range(d2_ints.shape[0]-1,1,-1):
      #print curr_scale
      rets[curr_scale]=[]
      for curr_x_in_scale in range(d2_ints.shape[1]):
         if d2_ints[curr_scale,curr_x_in_scale]!=0:
            #return trace_contour(curr_scale,curr_x_in_scale,d2_ints)
            x_intersect=trace_contour(curr_scale,curr_x_in_scale,d2_ints)
            if not x_intersect in x_intersect_to_scale:
               rets[curr_scale].append(x_intersect)
               x_intersect_to_scale[x_intersect] = curr_scale
            #print rets[curr_scale][-1]
    return rets,x_intersect_to_scale

"""
get means and medians is 

"""
def get_cp_medians(np.ndarray[np.float32_t,ndim=1] cps, np.ndarray[np.int64_t,ndim=1] edges):
   cdef np.ndarray medians = np.zeros(edges.shape[0]-1,dtype=np.float)
   cdef int i   

   for i in range(edges.shape[0]-1):
      medians[i] = np.median(cps[edges[i]:edges[i+1]])    

   return medians

def get_cp_means(np.ndarray[np.float32_t,ndim=1] cps, np.ndarray[np.int64_t,ndim=1] edges):
   cdef np.ndarray means = np.zeros(edges.shape[0]-1,dtype=np.float)
   cdef int i   

   for i in range(edges.shape[0]-1):
      means[i] = np.mean(cps[edges[i]:edges[i+1]])    

   return means

#def iterative_merge_edges(np.ndarray[np.float32_t,ndim=1] cps,np.ndarray[np.int64_t,ndim=1] edges, float min_merge_dif):
def iterative_merge_edges(np.ndarray[np.float32_t,ndim=1] cps,np.ndarray[np.int64_t,ndim=1] edges, float min_merge_dif):

   cdef float curr_min = 100.0
   cdef np.ndarray medians
   cdef np.ndarray deltas
   cdef np.ndarray merged_edges

   cdef int min_index


   while 1:
      medians=get_cp_medians(cps,np.array(edges))
      deltas=np.abs(np.diff(medians))   
      min_index=np.argmin(deltas)
      curr_min=deltas[min_index]
      if curr_min<min_merge_dif:
         merged_edges = np.zeros(edges.shape[0]-1,dtype=edges.dtype)
         merged_edges[:min_index+1]=edges[0:min_index+1]
         merged_edges[min_index+1:]=edges[min_index+2:]
         edges=merged_edges

      else:
         break

   return edges

      
      
   

   








