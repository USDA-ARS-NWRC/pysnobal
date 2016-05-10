"""
C implementation of Snobal

Run iSnobal for a single time step given an
array of Snobal instances, and dicts of inputs
"""


import cython
import numpy as np
cimport numpy as np
from cython.parallel import prange, parallel
from snobal import snobal
from libc.math cimport floor


@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def isnobal(np.ndarray[object, mode="c", ndim=2] s, 
               input1,
               input2,
               int nthreads=1):
    '''
    Call the function krige_grid in krige.c which will iterate over the grid
    within the C code
    
    Args:
        ta, tw, z, skvfac
    Out:
        thermal changed in place
    
    20160510 Scott Havens
    '''
    
    cdef int ny = s.shape[0]
    cdef int nx = s.shape[1]
    cdef int N = nx * ny
    cdef int index, ix, iy

    # convert the ta array to C
    cdef np.ndarray[object, mode="c", ndim=2] s_arr
    s_arr = np.ascontiguousarray(s, dtype=np.object)
    
#     with nogil, parallel():
    for index in range(N):
            
        # index = np.unravel_index(i, (ny,nx))
        if index < nx:
            ix = index
            iy = 0
        else:
            iy = int(floor(index / nx))
            ix = index - nx * iy
            
        s[iy,ix]

                
        
#         if s[index] is not None:
#         
#             in1 = {key: input1[key][index] for key in input1.keys()}
#             in2 = {key: input2[key][index] for key in input2.keys()}
#             
#             s[index].do_data_tstep(in1, in2)
    
    
    
    
#     cdef int ngrid
#     ngrid = ta.shape[0] * ta.shape[1]
#     
#     # convert the ta array to C
#     cdef np.ndarray[double, mode="c", ndim=2] ta_arr
#     ta_arr = np.ascontiguousarray(ta, dtype=np.float64)
#      
#     # convert the tw array to C
#     cdef np.ndarray[double, mode="c", ndim=2] tw_arr
#     tw_arr = np.ascontiguousarray(tw, dtype=np.float64)
#      
#     # convert the ta array to C
#     cdef np.ndarray[double, mode="c", ndim=2] z_arr
#     z_arr = np.ascontiguousarray(z, dtype=np.float64)
#      
#     # convert the skvfac to C
#     cdef np.ndarray[double, mode="c", ndim=2] skvfac_arr
#     skvfac_arr = np.ascontiguousarray(skvfac, dtype=np.float64)
#               
#     # call the C function
#     topotherm(ngrid, &ta_arr[0,0], &tw_arr[0,0], &z_arr[0,0], &skvfac_arr[0,0], nthreads, &thermal[0,0])

    return None
