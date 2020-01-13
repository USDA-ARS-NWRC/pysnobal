"""
Wrapper functions to the C function in libsnobal

20161010 Scott Havens
"""
import cython
from cython.parallel import prange, parallel
import numpy as np
cimport numpy as np

from libc.stdlib cimport abort, calloc, malloc
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF


# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

cdef extern from "time.h":
    ctypedef unsigned long clock_t
    cdef clock_t clock()
    cdef enum:
        CLOCKS_PER_SEC

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

    cdef TSTEP_REC tstep_info[4]

    int precip_now;
    double m_pp;
    double percent_snow;
    double rho_snow;
    double T_pp;

    int layer_count;
    double elevation;
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

# ctypedef struct OUTPUT_REC:
#         int masked;
#         double elevation;
#         double z_0;
#         double z_s;
#         double rho;
#         double T_s_0;
#         double T_s_l;
#         double T_s;
#         double h2o_sat;
#         int layer_count;
#         double R_n_bar;
#         double H_bar;
#         double L_v_E_bar;
#         double G_bar;
#         double G_0_bar;
#         double M_bar;
#         double delta_Q_bar;
#         double delta_Q_0_bar;
#         double E_s_sum;
#         double melt_sum;
#         double ro_pred_sum;

cdef extern from "pysnobal.h":
    #cdef int call_snobal(int N, int nthreads, int first_step, TSTEP_REC tstep_info[4], OUTPUT_REC** output_rec, INPUT_REC_ARR* input1, INPUT_REC_ARR* input2, PARAMS params, OUTPUT_REC_ARR* output1);
    cdef int call_snobal(int N, int nthreads, int first_step, TSTEP_REC tstep_info[4], INPUT_REC_ARR* input1, INPUT_REC_ARR* input2, PARAMS params, OUTPUT_REC_ARR* output1);

    ctypedef struct OUTPUT_REC:
        int masked;
        double current_time;
        double time_since_out;
        double elevation;
        double z_0;
        double rho;
        double T_s_0;
        double T_s_l;
        double T_s;
        double h2o_sat;
        double h2o_max;
        double h2o;
        double h2o_vol;
        double h2o_total;
        int layer_count;
        double cc_s_0;
        double cc_s_l;
        double cc_s;
        double m_s_0;
        double m_s_l;
        double m_s;
        double z_s_0;
        double z_s_l;
        double z_s;
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

    ctypedef struct OUTPUT_REC_ARR:
        int* masked;
        double* current_time;
        double* time_since_out;
        double* elevation;
        double* z_0;
        double* rho;
        double* T_s_0;
        double* T_s_l;
        double* T_s;
        double* h2o_sat;
        double* h2o_max;
        double* h2o;
        double* h2o_vol;
        double* h2o_total;
        int* layer_count;
        double* cc_s_0;
        double* cc_s_l;
        double* cc_s;
        double* m_s_0;
        double* m_s_l;
        double* m_s;
        double* z_s_0;
        double* z_s_l;
        double* z_s;
        double* R_n_bar;
        double* H_bar;
        double* L_v_E_bar;
        double* G_bar;
        double* G_0_bar;
        double* M_bar;
        double* delta_Q_bar;
        double* delta_Q_0_bar;
        double* E_s_sum;
        double* melt_sum;
        double* ro_pred_sum;

    ctypedef struct INPUT_REC_ARR:
        double* S_n;
        double* I_lw;
        double* T_a;
        double* e_a;
        double* u;
        double* T_g;
        double* m_pp;
        double* percent_snow;
        double* rho_snow;
        double* T_pp;

    ctypedef struct PARAMS:
        double z_u;
        double z_T;
        double z_g;
        int relative_heights;
        double max_h2o_vol;
        double max_z_s_0;



@cython.boundscheck(False)
@cython.wraparound(False)
# https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
def do_tstep_grid(input1, input2, output_rec, tstep_rec, mh, params, int first_step=1, int nthreads=1):
    """
    Do the timestep given the inputs, model state, and measurement heights
    There is no first_step value since the snow state records were already
    pulled in an initialized. Therefore only the values need to be pulled
    out before calling 'init_snow()'
    """
    #cdef int N = len(output_rec['elevation'])
    cdef int N = (output_rec['elevation']).size
    cdef int n
    shp = output_rec['elevation'].shape
    # timing stuff
    cdef clock_t start1, end1, start2, end2, start3, end3
    cdef double cpu_time_used1, cpu_time_used2, cpu_time_used3

    # measurement heights and parameters
    cdef PARAMS c_params
    c_params.z_u = mh['z_u']
    c_params.z_T = mh['z_t']
    c_params.z_g = mh['z_g']
    c_params.relative_heights = int(params['relative_heights'])
    c_params.max_h2o_vol = params['max_h2o_vol']
    c_params.max_z_s_0 = params['max_z_s_0']
#


#     cdef TSTEP_REC tstep_info[4]
#     cdef TSTEP_REC **tstep_info = <TSTEP_REC**> PyMem_Malloc(4 * sizeof(TSTEP_REC*));
#     for n in range(N):
#         tstep_info[n] = <TSTEP_REC *>PyMem_Malloc(sizeof(TSTEP_REC))
#         tstep_info[n].level = n

    for i in range(len(tstep_rec)):
        tstep_info[i].level = int(tstep_rec[i]['level'])
        if tstep_rec[i]['time_step'] is not None:
            tstep_info[i].time_step = tstep_rec[i]['time_step']
        if tstep_rec[i]['intervals'] is not None:
            tstep_info[i].intervals = int(tstep_rec[i]['intervals'])
        if tstep_rec[i]['threshold'] is not None:
            tstep_info[i].threshold = tstep_rec[i]['threshold']
        tstep_info[i].output = int(tstep_rec[i]['output'])

    # start1 = clock()
    cdef OUTPUT_REC_ARR output1_c

    cdef np.ndarray[double, mode="c", ndim=2] output1_z_0
    output1_z_0 = np.ascontiguousarray(output_rec['z_0'], dtype=np.float64)
    output1_c.z_0 = &output1_z_0[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_z_s_0
    output1_z_s_0 = np.ascontiguousarray(output_rec['z_s_0'], dtype=np.float64)
    output1_c.z_s_0 = &output1_z_s_0[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_current_time
    output1_current_time = np.ascontiguousarray(output_rec['current_time'], dtype=np.float64)
    output1_c.current_time = &output1_current_time[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_time_since_out
    output1_time_since_out = np.ascontiguousarray(output_rec['time_since_out'], dtype=np.float64)
    output1_c.time_since_out = &output1_time_since_out[0,0]

    cdef np.ndarray[int, mode="c", ndim=2] output1_masked
    output1_masked = np.ascontiguousarray(output_rec['mask'], dtype=np.int32)
    output1_c.masked = &output1_masked[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_elevation
    output1_elevation = np.ascontiguousarray(output_rec['elevation'], dtype=np.float64)
    output1_c.elevation = &output1_elevation[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_z_s_l
    output1_z_s_l = np.ascontiguousarray(output_rec['z_s_l'], dtype=np.float64)
    output1_c.z_s_l = &output1_z_s_l[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_z_s
    output1_z_s = np.ascontiguousarray(output_rec['z_s'], dtype=np.float64)
    output1_c.z_s = &output1_z_s[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_rho
    output1_rho = np.ascontiguousarray(output_rec['rho'], dtype=np.float64)
    output1_c.rho = &output1_rho[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_T_s_0
    output1_T_s_0 = np.ascontiguousarray(output_rec['T_s_0'], dtype=np.float64)
    output1_c.T_s_0 = &output1_T_s_0[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_T_s_l
    output1_T_s_l = np.ascontiguousarray(output_rec['T_s_l'], dtype=np.float64)
    output1_c.T_s_l = &output1_T_s_l[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_T_s
    output1_T_s = np.ascontiguousarray(output_rec['T_s'], dtype=np.float64)
    output1_c.T_s = &output1_T_s[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_cc_s_0
    output1_cc_s_0 = np.ascontiguousarray(output_rec['cc_s_0'], dtype=np.float64)
    output1_c.cc_s_0 = &output1_cc_s_0[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_cc_s_l
    output1_cc_s_l = np.ascontiguousarray(output_rec['cc_s_l'], dtype=np.float64)
    output1_c.cc_s_l = &output1_cc_s_l[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_cc_s
    output1_cc_s = np.ascontiguousarray(output_rec['cc_s'], dtype=np.float64)
    output1_c.cc_s = &output1_cc_s[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_m_s_0
    output1_m_s_0 = np.ascontiguousarray(output_rec['m_s_0'], dtype=np.float64)
    output1_c.m_s_0 = &output1_m_s_0[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_m_s_l
    output1_m_s_l = np.ascontiguousarray(output_rec['m_s_l'], dtype=np.float64)
    output1_c.m_s_l = &output1_m_s_l[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_m_s
    output1_m_s = np.ascontiguousarray(output_rec['m_s'], dtype=np.float64)
    output1_c.m_s = &output1_m_s[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_h2o_sat
    output1_h2o_sat = np.ascontiguousarray(output_rec['h2o_sat'], dtype=np.float64)
    output1_c.h2o_sat = &output1_h2o_sat[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_h2o_max
    output1_h2o_max = np.ascontiguousarray(output_rec['h2o_max'], dtype=np.float64)
    output1_c.h2o_max = &output1_h2o_max[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_h2o
    output1_h2o = np.ascontiguousarray(output_rec['h2o'], dtype=np.float64)
    output1_c.h2o = &output1_h2o[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_h2o_vol
    output1_h2o_vol = np.ascontiguousarray(output_rec['h2o_vol'], dtype=np.float64)
    output1_c.h2o_vol = &output1_h2o_vol[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_h2o_total
    output1_h2o_total = np.ascontiguousarray(output_rec['h2o_total'], dtype=np.float64)
    output1_c.h2o_total = &output1_h2o_total[0,0]

    cdef np.ndarray[int, mode="c", ndim=2] output1_layer_count
    output1_layer_count = np.ascontiguousarray(output_rec['layer_count'], dtype=np.int32)
    output1_c.layer_count = &output1_layer_count[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_R_n_bar
    output1_R_n_bar = np.ascontiguousarray(output_rec['R_n_bar'], dtype=np.float64)
    output1_c.R_n_bar = &output1_R_n_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_H_bar
    output1_H_bar = np.ascontiguousarray(output_rec['H_bar'], dtype=np.float64)
    output1_c.H_bar = &output1_H_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_L_v_E_bar
    output1_L_v_E_bar = np.ascontiguousarray(output_rec['L_v_E_bar'], dtype=np.float64)
    output1_c.L_v_E_bar = &output1_L_v_E_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_G_bar
    output1_G_bar = np.ascontiguousarray(output_rec['G_bar'], dtype=np.float64)
    output1_c.G_bar = &output1_G_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_G_0_bar
    output1_G_0_bar = np.ascontiguousarray(output_rec['G_0_bar'], dtype=np.float64)
    output1_c.G_0_bar = &output1_G_0_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_M_bar
    output1_M_bar = np.ascontiguousarray(output_rec['M_bar'], dtype=np.float64)
    output1_c.M_bar = &output1_M_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_delta_Q_bar
    output1_delta_Q_bar = np.ascontiguousarray(output_rec['delta_Q_bar'], dtype=np.float64)
    output1_c.delta_Q_bar = &output1_delta_Q_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_delta_Q_0_bar
    output1_delta_Q_0_bar = np.ascontiguousarray(output_rec['delta_Q_0_bar'], dtype=np.float64)
    output1_c.delta_Q_0_bar = &output1_delta_Q_0_bar[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_E_s_sum
    output1_E_s_sum = np.ascontiguousarray(output_rec['E_s_sum'], dtype=np.float64)
    output1_c.E_s_sum = &output1_E_s_sum[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_melt_sum
    output1_melt_sum = np.ascontiguousarray(output_rec['melt_sum'], dtype=np.float64)
    output1_c.melt_sum = &output1_melt_sum[0,0]

    cdef np.ndarray[double, mode="c", ndim=2] output1_ro_pred_sum
    output1_ro_pred_sum = np.ascontiguousarray(output_rec['ro_pred_sum'], dtype=np.float64)
    output1_c.ro_pred_sum = &output1_ro_pred_sum[0,0]

    # end1 = clock()
    # cpu_time_used1 = (<double> (end1 - start1)) / CLOCKS_PER_SEC
    # print('time 1 {}'.format(cpu_time_used1))
    # start15 = clock()
    #------------------------------------------------------------------------------
    # PREPARE INPUT1 FOR C

    cdef INPUT_REC_ARR input1_c

    # convert the S_n to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_Sn
    input1_Sn = np.ascontiguousarray(input1['S_n'], dtype=np.float64)
    input1_c.S_n = &input1_Sn[0,0]

    # convert the I_lw to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_I_lw
    input1_I_lw = np.ascontiguousarray(input1['I_lw'], dtype=np.float64)
    input1_c.I_lw = &input1_I_lw[0,0]

    # convert the T_a to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_Ta
    input1_Ta = np.ascontiguousarray(input1['T_a'], dtype=np.float64)
    input1_c.T_a = &input1_Ta[0,0] # For some reason this isn't needed, most likely b/c numpy has already allocated it    input1_c.T_a = <double *> PyMem_Malloc(N * sizeof(double))

    # convert the e_a to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_e_a
    input1_e_a = np.ascontiguousarray(input1['e_a'], dtype=np.float64)
    input1_c.e_a = &input1_e_a[0,0]

    # convert the u to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_u
    input1_u = np.ascontiguousarray(input1['u'], dtype=np.float64)
    input1_c.u = &input1_u[0,0]

    # convert the T_g to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_T_g
    input1_T_g = np.ascontiguousarray(input1['T_g'], dtype=np.float64)
    input1_c.T_g = &input1_T_g[0,0]

    # convert the m_pp to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_m_pp
    input1_m_pp = np.ascontiguousarray(input1['m_pp'], dtype=np.float64)
    input1_c.m_pp = &input1_m_pp[0,0]

    # convert the percent_snow to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_percent_snow
    input1_percent_snow = np.ascontiguousarray(input1['percent_snow'], dtype=np.float64)
    input1_c.percent_snow = &input1_percent_snow[0,0]

    # convert the rho_snow to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_rho_snow
    input1_rho_snow = np.ascontiguousarray(input1['rho_snow'], dtype=np.float64)
    input1_c.rho_snow = &input1_rho_snow[0,0]

    # convert the T_pp to C
    cdef np.ndarray[double, mode="c", ndim=2] input1_T_pp
    input1_T_pp = np.ascontiguousarray(input1['T_pp'], dtype=np.float64)
    input1_c.T_pp = &input1_T_pp[0,0]


    #------------------------------------------------------------------------------
    # PREPARE INPUT2 FOR C

    cdef INPUT_REC_ARR input2_c

    # convert the S_n to C
    cdef np.ndarray[double, mode="c", ndim=2] input2_Sn
    input2_Sn = np.ascontiguousarray(input2['S_n'], dtype=np.float64)
    input2_c.S_n = &input2_Sn[0,0]

    # convert the I_lw to C
    cdef np.ndarray[double, mode="c", ndim=2] input2_I_lw
    input2_I_lw = np.ascontiguousarray(input2['I_lw'], dtype=np.float64)
    input2_c.I_lw = &input2_I_lw[0,0]

    # convert the T_a to C
    cdef np.ndarray[double, mode="c", ndim=2] input2_Ta
    input2_Ta = np.ascontiguousarray(input2['T_a'], dtype=np.float64)
    input2_c.T_a = &input2_Ta[0,0]

    # convert the e_a to C
    cdef np.ndarray[double, mode="c", ndim=2] input2_e_a
    input2_e_a = np.ascontiguousarray(input2['e_a'], dtype=np.float64)
    input2_c.e_a = &input2_e_a[0,0]

    # convert the u to C
    cdef np.ndarray[double, mode="c", ndim=2] input2_u
    input2_u = np.ascontiguousarray(input2['u'], dtype=np.float64)
    input2_c.u = &input2_u[0,0]

    # convert the T_g to C
    cdef np.ndarray[double, mode="c", ndim=2] input2_T_g
    input2_T_g = np.ascontiguousarray(input2['T_g'], dtype=np.float64)
    input2_c.T_g = &input2_T_g[0,0]

    # timing stuff
    # end15 = clock()
    # cpu_time_used15 = (<double> (end15 - start15)) / CLOCKS_PER_SEC
    # print('time 1.5 {}'.format(cpu_time_used15))
    # start2 = clock()

    #------------------------------------------------------------------------------
    # Call the model
    # rt = call_snobal(N, nthreads, first_step, tstep_info, out_c, &input1_c, &input2_c, c_params, &output1_c)
    rt = call_snobal(N, nthreads, first_step, tstep_info, &input1_c, &input2_c, c_params, &output1_c)
    if rt != -1:
        return rt

    # end2 = clock()
    # cpu_time_used2 = (<double> (end2 - start2)) / CLOCKS_PER_SEC
    # print('time 2 {}'.format(cpu_time_used2))
    # start3 = clock()

    cdef np.npy_intp shp_np[2]
    shp_np[:] = (shp[0], shp[1])
    output_rec['z_0'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.z_0)
    output_rec['z_s_0'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.z_s_0)

    output_rec['current_time'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.current_time)
    output_rec['time_since_out'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.time_since_out)

    output_rec['elevation'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.elevation)
    output_rec['rho'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.rho)
    output_rec['T_s_0'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.T_s_0)
    output_rec['T_s_l'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.T_s_l)
    output_rec['T_s'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.T_s)
    output_rec['h2o_sat'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.h2o_sat)
    output_rec['h2o_max'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.h2o_max)
    output_rec['h2o'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.h2o)
    output_rec['h2o_vol'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.h2o_vol)
    output_rec['h2o_total'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.h2o_total)
    output_rec['layer_count'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_INT32, output1_c.layer_count)
    output_rec['cc_s_0'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.cc_s_0)
    output_rec['cc_s_l'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.cc_s_l)
    output_rec['cc_s'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.cc_s)
    output_rec['m_s_0'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.m_s_0)
    output_rec['m_s_l'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.m_s_l)
    output_rec['m_s'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.m_s)
    output_rec['z_s_l'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.z_s_l)
    output_rec['z_s'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.z_s)

    output_rec['R_n_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.R_n_bar)
    output_rec['H_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.H_bar)
    output_rec['L_v_E_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.L_v_E_bar)
    output_rec['G_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.G_bar)
    output_rec['G_0_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.G_0_bar)
    output_rec['M_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.M_bar)
    output_rec['delta_Q_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.delta_Q_bar)
    output_rec['delta_Q_0_bar'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.delta_Q_0_bar)
    output_rec['E_s_sum'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.E_s_sum)
    output_rec['melt_sum'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.melt_sum)
    output_rec['ro_pred_sum'][:] = np.PyArray_SimpleNewFromData(2, shp_np, np.NPY_FLOAT64, output1_c.ro_pred_sum)

    # end3 = clock()
    # cpu_time_used3 = (<double> (end3 - start3)) / CLOCKS_PER_SEC
    # print('time 3 {}'.format(cpu_time_used3))

    return rt




# We need to build an array-wrapper class to deallocate our array when
# the Python object is deleted.
# From https://gist.github.com/GaelVaroquaux/1249305
cdef class ArrayWrapper:
    cdef void* data_ptr
    cdef int size

    cdef set_data(self, int size, void* data_ptr):
        """ Set the data of the array
        This cannot be done in the constructor as it must recieve C-level
        arguments.
        Parameters:
        -----------
        size: int
            Length of the array.
        data_ptr: void*
            Pointer to the data
        """
        self.data_ptr = data_ptr
        self.size = size

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef np.npy_intp shape[1]
        shape[0] = <np.npy_intp> self.size
        # Create a 1D array, of length 'size'
        ndarray = np.PyArray_SimpleNewFromData(1, shape,
                                               np.NPY_INT, self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def do_tstep(input1, input2, output_rec, tstep_rec, mh, params, first_step=True):
    """
    Do the timestep given the inputs, model state, and measurement heights
    There is no first_step value since the snow state records were already
    pulled in an initialized. Therefore only the values need to be pulled
    out before calling 'init_snow()'
    """

    cdef int N = len(output_rec['elevation'])
#     cdef TSTEP_REC tstep_info[4]

    for i in range(len(tstep_rec)):
        tstep_info[i].level = int(tstep_rec[i]['level'])
        if tstep_rec[i]['time_step'] is not None:
            tstep_info[i].time_step = tstep_rec[i]['time_step']
        if tstep_rec[i]['intervals'] is not None:
            tstep_info[i].intervals = int(tstep_rec[i]['intervals'])
        if tstep_rec[i]['threshold'] is not None:
            tstep_info[i].threshold = tstep_rec[i]['threshold']
        tstep_info[i].output = int(tstep_rec[i]['output'])



    # array of pointers to OUTPUT_REC
#     cdef OUTPUT_REC **out_c = <OUTPUT_REC**> malloc(N * sizeof(OUTPUT_REC*));
#     cdef int n
#     cdef OUTPUT_REC tmp
#     for n in range(N):
# #         tmp =
#         out_c[n] = malloc(sizeof(OUTPUT_REC))    # initialize the memory (in heap) at that pointer




    # loop through the grid
    rt = True
    for (i,j), z in np.ndenumerate(output_rec['elevation']):
#     with gil, parallel(num_threads=4):
#         for n in prange(N):

        # extract_data.c
        #check to see if point is masked, since 1=run point, it's "not" masked
        masked = output_rec['mask'][i,j]
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
            current_time = output_rec['current_time'][i,j]
            time_since_out = output_rec['time_since_out'][i,j]

            # measurement heights and parameters
            z_u = mh['z_u']
            z_T = mh['z_t']
            z_g = mh['z_g']
            relative_heights = int(params['relative_heights'])
            max_h2o_vol = params['max_h2o_vol']
            max_z_s_0 = params['max_z_s_0']

            # get the input records
            input_rec1.I_lw = input1['I_lw'][i,j]
            input_rec1.T_a  = input1['T_a'][i,j]
            input_rec1.e_a  = input1['e_a'][i,j]
            input_rec1.u    = input1['u'][i,j]
            input_rec1.T_g  = input1['T_g'][i,j]
            input_rec1.S_n  = input1['S_n'][i,j]

            input_rec2.I_lw = input2['I_lw'][i,j]
            input_rec2.T_a  = input2['T_a'][i,j]
            input_rec2.e_a  = input2['e_a'][i,j]
            input_rec2.u    = input2['u'][i,j]
            input_rec2.T_g  = input2['T_g'][i,j]
            input_rec2.S_n  = input2['S_n'][i,j]

            m_pp         = input1['m_pp'][i,j]
            percent_snow = input1['percent_snow'][i,j]
            rho_snow     = input1['rho_snow'][i,j]
            T_pp         = input1['T_pp'][i,j]

            precip_now = 0
            if m_pp > 0:
                precip_now = 1

            # get the model state
            elevation       = output_rec['elevation'][i,j]
            z_0             = output_rec['z_0'][i,j]
            z_s             = output_rec['z_s'][i,j]
            rho             = output_rec['rho'][i,j]
            T_s_0           = output_rec['T_s_0'][i,j]
            T_s_l           = output_rec['T_s_l'][i,j]
            T_s             = output_rec['T_s'][i,j]
            h2o_sat         = output_rec['h2o_sat'][i,j]
            layer_count     = output_rec['layer_count'][i,j]

            R_n_bar         = output_rec['R_n_bar'][i,j]
            H_bar           = output_rec['H_bar'][i,j]
            L_v_E_bar       = output_rec['L_v_E_bar'][i,j]
            G_bar           = output_rec['G_bar'][i,j]
            G_0_bar         = output_rec['G_0_bar'][i,j]
            M_bar           = output_rec['M_bar'][i,j]
            delta_Q_bar     = output_rec['delta_Q_bar'][i,j]
            delta_Q_0_bar   = output_rec['delta_Q_0_bar'][i,j]
            E_s_sum         = output_rec['E_s_sum'][i,j]
            melt_sum        = output_rec['melt_sum'][i,j]
            ro_pred_sum     = output_rec['ro_pred_sum'][i,j]

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

            # do_data_tstep.c
            dt = do_data_tstep()
            if dt == 0:
                rt = False
                if N > 1:
                    abort()
                else:
                    break

            output_rec['current_time'][i,j] = current_time
            output_rec['time_since_out'][i,j] = time_since_out

            output_rec['elevation'][i,j] = elevation
            output_rec['z_0'][i,j] = z_0
            output_rec['rho'][i,j] = rho
            output_rec['T_s_0'][i,j] = T_s_0
            output_rec['T_s_l'][i,j] = T_s_l
            output_rec['T_s'][i,j] = T_s
            output_rec['h2o_sat'][i,j] = h2o_sat
            output_rec['h2o_max'][i,j] = h2o_max
            output_rec['h2o'][i,j] = h2o
            output_rec['layer_count'][i,j] = layer_count
            output_rec['cc_s_0'][i,j] = cc_s_0
            output_rec['cc_s_l'][i,j] = cc_s_l
            output_rec['cc_s'][i,j] = cc_s
            output_rec['m_s_0'][i,j] = m_s_0
            output_rec['m_s_l'][i,j] = m_s_l
            output_rec['m_s'][i,j] = m_s
            output_rec['z_s_0'][i,j] = z_s_0
            output_rec['z_s_l'][i,j] = z_s_l
            output_rec['z_s'][i,j] = z_s

            output_rec['R_n_bar'][i,j] = R_n_bar
            output_rec['H_bar'][i,j] = H_bar
            output_rec['L_v_E_bar'][i,j] = L_v_E_bar
            output_rec['G_bar'][i,j] = G_bar
            output_rec['G_0_bar'][i,j] = G_0_bar
            output_rec['M_bar'][i,j] = M_bar
            output_rec['delta_Q_bar'][i,j] = delta_Q_bar
            output_rec['delta_Q_0_bar'][i,j] = delta_Q_0_bar
            output_rec['E_s_sum'][i,j] = E_s_sum
            output_rec['melt_sum'][i,j] = melt_sum
            output_rec['ro_pred_sum'][i,j] = ro_pred_sum

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