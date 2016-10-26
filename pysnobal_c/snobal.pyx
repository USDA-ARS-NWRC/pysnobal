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
    
    ctypedef struct INPUT_REC:
        double S_n;
        double I_lw;
        double T_a;
        double e_a;
        double u;
        double T_g;
        double ro;
    
    INPUT_REC  input_rec1;
    INPUT_REC  input_rec2;
    
    ctypedef struct TSTEP_REC:
        int level;
        double time_step;
        int intervals;
        double threshold;
        int output;
    
    TSTEP_REC tstep_info[4];
    
    int precip_now;
    double m_pp;
    double percent_snow;
    double rho_snow;
    double T_pp;
    
    int layer_count;
    double z_0;
    double rho;
    double T_s;
    double T_s_0;
    double T_s_l;
    double h2o_sat;
    double h2o;
    double h2o_max;
    double P_a;
    double m_s;
    double m_s_0;
    double m_s_l;
    double cc_s;
    double cc_s_0;
    double cc_s_l;
    double z_s;
    double z_s_0;
    double z_s_l;
        
    double R_n_bar;
    double H_bar;
    double L_v_E_bar;
    double G_bar;
    double G_0_bar;
    double M_bar;
    double delta_Q_bar;
    double delta_Q_0_bar;
    double E_s_sum;
    double melt_sum;
    double ro_pred_sum;
    
    double current_time;
    double time_since_out;
    
    int relative_hts;
    double z_g;
    double z_u;
    double z_T;
    double max_h2o_vol;
    double max_z_s_0;


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

def do_tstep(input1, input2, output_rec, tstep_rec, mh, params, first_step=False):
    """
    Do the timestep given the inputs, model state, and measurement heights
    There is no first_step value since the snow state records were already
    pulled in an initialized. Therefore only the values need to be pulled 
    out before calling 'init_snow()'
    """
    
#     n = 0
#     n = np.unravel_index(n, output_rec['elevation'].shape)
    
    for i in range(len(tstep_rec)):
        tstep_info[i].level = int(tstep_rec[i]['level'])
        if tstep_rec[i]['time_step'] is not None:
            tstep_info[i].time_step = tstep_rec[i]['time_step']
        if tstep_rec[i]['intervals'] is not None:
            tstep_info[i].intervals = int(tstep_rec[i]['intervals'])
        if tstep_rec[i]['threshold'] is not None:
            tstep_info[i].threshold = tstep_rec[i]['threshold']
        tstep_info[i].output = int(tstep_rec[i]['output'])
    
    # loop through the grid
    for n, z in np.ndenumerate(output_rec['elevation']):
            
        # extract_data.c
        #check to see if point is masked, since 1=run point, it's "not" masked
        masked = output_rec['mask'][n]
        if masked:
            
            # since snobal use global variables extensively, here is the ugly
            # interface with Python
            global current_time, time_since_out
            global m_pp, percent_snow, rho_snow, T_pp, precip_now
            global z_0, rho, T_s_0, T_s_l, T_s, h2o_sat, h2o_max, layer_count, P_a, h2o
            global m_s_0, m_s_l, m_s, cc_s_0, cc_s_l, cc_s, z_s_0, z_s_l, z_s
            global z_u, z_T, z_g, relative_heights, max_h2o_vol, max_z_s_0
            global R_n_bar, H_bar, L_v_E_bar, G_bar, G_0_bar, M_bar, delta_Q_bar, delta_Q_0_bar
            global E_s_sum, melt_sum, ro_pred_sum 
             
            # time variables
            current_time = output_rec['current_time'][n]
            time_since_out = output_rec['time_since_out'][n]
            
            # measurement heights and parameters
            z_u = mh['z_u']
            z_T = mh['z_t']
            z_g = mh['z_g']
            relative_heights = int(params['relative_heights'])
            max_h2o_vol = params['max_h2o_vol']
            max_z_s_0 = params['max_z_s_0']
             
            # get the input records
            input_rec1.I_lw = input1['I_lw'][n]
            input_rec1.T_a  = input1['T_a'][n]
            input_rec1.e_a  = input1['e_a'][n]
            input_rec1.u    = input1['u'][n]
            input_rec1.T_g  = input1['T_g'][n] 
            input_rec1.S_n  = input1['S_n'][n]
             
            input_rec2.I_lw = input2['I_lw'][n]
            input_rec2.T_a  = input2['T_a'][n]
            input_rec2.e_a  = input2['e_a'][n]
            input_rec2.u    = input2['u'][n]
            input_rec2.T_g  = input2['T_g'][n] 
            input_rec2.S_n  = input2['S_n'][n]
            
            m_pp         = input1['m_pp'][n]
            percent_snow = input1['percent_snow'][n]
            rho_snow     = input1['rho_snow'][n]
            T_pp         = input1['T_pp'][n]
            
            precip_now = 0
            if m_pp > 0:
                precip_now = 1
             
            # get the model state
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
            G_0_bar         = output_rec['G_0_bar'][n]
            M_bar           = output_rec['M_bar'][n]
            delta_Q_bar     = output_rec['delta_Q_bar'][n]
            delta_Q_0_bar   = output_rec['delta_Q_0_bar'][n]
            E_s_sum         = output_rec['E_s_sum'][n]
            melt_sum        = output_rec['melt_sum'][n]
            ro_pred_sum     = output_rec['ro_pred_sum'][n]
     
    #         print z_0
     
            # establish conditions for snowpack
            # the firs step mimic's snobal which only calls init_snow once. This
            # might mean that the model states will need to be saved in memory
            # or there will be a slight descrepancy with isnobal. But with this,
            # there should be a descrepancy in isnobal as well
            if first_step:
                init_snow()
     
            # set air pressure from site elev
            P_a = HYSTAT(SEA_LEVEL, STD_AIRTMP, STD_LAPSE, (elevation / 1000.0),
                GRAVITY, MOL_AIR)
         
    #         print cc_s_0
         
            # do_data_tstep.c
            do_data_tstep()
         
         
            # assign_buffers.c
    #         print R_n_bar
        
            output_rec['current_time'][n] = current_time
            output_rec['time_since_out'][n] = time_since_out
    
            output_rec['elevation'][n] = elevation
            output_rec['z_0'][n] = z_0
            output_rec['rho'][n] = rho
            output_rec['T_s_0'][n] = T_s_0
            output_rec['T_s_l'][n] = T_s_l
            output_rec['T_s'][n] = T_s
            output_rec['h2o_sat'][n] = h2o_sat
            output_rec['h2o_max'][n] = h2o_max
            output_rec['h2o'][n] = h2o
            output_rec['layer_count'][n] = layer_count
            output_rec['cc_s_0'][n] = cc_s_0
            output_rec['cc_s_l'][n] = cc_s_l
            output_rec['cc_s'][n] = cc_s
            output_rec['m_s_0'][n] = m_s_0
            output_rec['m_s_l'][n] = m_s_l
            output_rec['m_s'][n] = m_s
            output_rec['z_s_0'][n] = z_s_0
            output_rec['z_s_l'][n] = z_s_l
            output_rec['z_s'][n] = z_s
     
            output_rec['R_n_bar'][n] = R_n_bar
            output_rec['H_bar'][n] = H_bar
            output_rec['L_v_E_bar'][n] = L_v_E_bar
            output_rec['G_bar'][n] = G_bar
            output_rec['G_0_bar'][n] = G_0_bar
            output_rec['M_bar'][n] = M_bar
            output_rec['delta_Q_bar'][n] = delta_Q_bar
            output_rec['delta_Q_0_bar'][n] = delta_Q_0_bar
            output_rec['E_s_sum'][n] = E_s_sum
            output_rec['melt_sum'][n] = melt_sum
            output_rec['ro_pred_sum'][n] = ro_pred_sum

    return None


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