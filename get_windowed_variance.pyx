import numpy as np
cimport numpy as np

#@cython.boundscheck(False)
#@cython.wraparound(False)

cdef int calc_var( np.ndarray[ np.double_t, ndim=1 ] cps , np.ndarray[ np.double_t, ndim=1 ] means, np.ndarray[ np.double_t, ndim=1 ] variances,  int half_width ): 
   
    cdef full_width = (2*half_width) + 1
    cdef int i,j,k
    cdef float squared_sum
    
    for i in range(cps.shape[0]):
        
        squared_sum = 0
        for k in range( -half_width, half_width+1, 1):
            j=i+k
            if j <0: j=0
            if j>= cps.shape[0]: j = cps.shape[0]-1
            squared_sum += (cps[j]-means[i])*(cps[j]-means[i])
        variances[i] = squared_sum / full_width

def get_windowed_variance( np.ndarray[ np.double_t, ndim=1 ] cps , int half_width ):

    cdef np.ndarray variances = np.zeros(cps.shape[0], dtype = np.float)
    cdef np.ndarray means = np.zeros(cps.shape[0], dtype = np.float)
    
    cdef full_width = (2*half_width) + 1
    cdef int ret 
    means = np.convolve( cps, np.ones( full_width ) ) / full_width
    
    ret =  calc_var(cps, means, variances, half_width)
    
    return variances
    

   
