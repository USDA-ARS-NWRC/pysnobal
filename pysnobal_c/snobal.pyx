"""
Wrapper functions to the C function in libsnobal

20161010 Scott Havens
"""

import cython
import numpy as np
cimport numpy as np

# from libc.stdlib cimport free
# from cpython cimport PyObject, Py_INCREF

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()


cdef extern from "snobal.h":
    void init_snow();
    int do_data_tstep();

cdef extern from "envphys.h":
    cdef double SEA_LEVEL;
    cdef double STD_AIRTMP;
    cdef double STD_LAPSE;
    cdef double GRAVITY;
    cdef double MOL_AIR;
    cdef double HYSTAT(double pb, double tb, double L, double h, double g, double m); 



def initialize(params, tstep_info, sn, mh):
    """
    Initialize the Snobal model given the input records
    
    Args:
        params:
        tstep_info:
        sn:
        mh:
    """
    
    return None

def do_tstep(input1, input2, output_rec, mh):
    """
    Do the timestep given the inputs, model state, and measurement heights
    There is no first_step value since the snow state records were already
    pulled in an initialized. Therefore only the values need to be pulled 
    out before calling 'init_snow()'
    """
    
    n = 0
    
    cdef double cc_s_0
    global cc_s_0
    
    # extract_data.c
    #check to see if point is masked
    masked = output_rec['masked'][n]
    if masked:
        rt = False

    else:
        elevation       = output_rec['elevation'][n]
        z_0             = output_rec['z_0'][n]
        z_s             = output_rec['z_s'][n]
        rho             = output_rec['rho'][n]
        T_s_0           = output_rec['T_s_0'][n]
        T_s_l           = output_rec['T_s_l'][n]
        T_s             = output_rec['T_s'][n]
        h2o_sat         = output_rec['h2o_sat'][n]
        layer_count     = output_rec['layer_count'][n]
 
        R_n_bar         = output_rec['R_n_bar'][n]
        H_bar           = output_rec['H_bar'][n]
        L_v_E_bar       = output_rec['L_v_E_bar'][n]
        G_bar           = output_rec['G_bar'][n]
        M_bar           = output_rec['M_bar'][n]
        delta_Q_bar     = output_rec['delta_Q_bar'][n]
        E_s_sum         = output_rec['E_s_sum'][n]
        melt_sum        = output_rec['melt_sum'][n]
        ro_pred_sum     = output_rec['ro_pred_sum'][n]

        # establish conditions for snowpack
        init_snow()

        # set air pressure from site elev
        P_a = HYSTAT(SEA_LEVEL, STD_AIRTMP, STD_LAPSE, (elevation / 1000.0),
            GRAVITY, MOL_AIR)
    
#         print cc_s_0
    
        # do_data_tstep.c
        do_data_tstep()
    
    
        # assign_buffers.c
        print R_n_bar
        output_rec['elevation'][n] = elevation
        output_rec['z_0'][n] = z_0
        output_rec['z_s'][n] = z_s
        output_rec['rho'][n] = rho
        output_rec['T_s_0'][n] = T_s_0
        output_rec['T_s_l'][n] = T_s_l
        output_rec['T_s'][n] = T_s
        output_rec['h2o_sat'][n] = h2o_sat
        output_rec['layer_count'][n] = layer_count

        output_rec['R_n_bar'][n] = R_n_bar
        output_rec['H_bar'][n] = H_bar
        output_rec['L_v_E_bar'][n] = L_v_E_bar
        output_rec['G_bar'][n] = G_bar
        output_rec['M_bar'][n] = M_bar
        output_rec['delta_Q_bar'][n] = delta_Q_bar
        output_rec['E_s_sum'][n] = E_s_sum
        output_rec['melt_sum'][n] = melt_sum
        output_rec['ro_pred_sum'][n] = ro_pred_sum

    return rt


@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def call_grid():
    '''
    Call the function krige_grid in krige.c which will iterate over the grid
    within the C code
    
    Args:
        ad - [nsta x nsta] matrix of distances between stations
        dgrid - [ngrid x nsta] matrix of distances between grid points and stations
        elevations - [nsta] array of station elevations
        weights (return) - [ngrid x nsta] matrix of kriging weights calculated
        nthreads - number of threads to use in parallel processing
        
    Out:
        weights changed in place
    
    20160222 Scott Havens
    '''
    
#     cdef int nsta, ngrid
#     ngrid, nsta = dgrid.shape[0], dgrid.shape[1]
#     
#     # convert the ad array to C
#     cdef np.ndarray[double, mode="c", ndim=2] ad_arr
#     ad_arr = np.ascontiguousarray(ad, dtype=np.float64)
#     
#     # convert the dgrid to C
#     cdef np.ndarray[double, mode="c", ndim=2] grid
#     grid = np.ascontiguousarray(dgrid, dtype=np.float64)
#             
#     # call the C function
#     krige_grid(nsta, ngrid, &ad_arr[0,0], &grid[0,0], &elevations[0], nthreads, &weights[0,0])

    return None