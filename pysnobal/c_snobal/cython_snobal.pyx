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


cdef extern from "pysnobal.h":

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
def do_tstep_point(input1, input2, output_rec, tstep_rec, mh, params, first_step=True):
    """
    Do the timestep given the inputs, model state, and measurement heights
    There is no first_step value since the snow state records were already
    pulled in an initialized. Therefore only the values need to be pulled
    out before calling 'init_snow()'
    """

    # cdef int N = len(output_rec['elevation'])
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
        
    rt = True

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
    current_time = output_rec['current_time']
    time_since_out = output_rec['time_since_out']

    # measurement heights and parameters
    z_u = mh['z_u']
    z_T = mh['z_t']
    z_g = mh['z_g']
    relative_heights = int(params['relative_heights'])
    max_h2o_vol = params['max_h2o_vol']
    max_z_s_0 = params['max_z_s_0']

    # get the input records
    input_rec1.I_lw = input1['I_lw']
    input_rec1.T_a  = input1['T_a']
    input_rec1.e_a  = input1['e_a']
    input_rec1.u    = input1['u']
    input_rec1.T_g  = input1['T_g']
    input_rec1.S_n  = input1['S_n']

    input_rec2.I_lw = input2['I_lw']
    input_rec2.T_a  = input2['T_a']
    input_rec2.e_a  = input2['e_a']
    input_rec2.u    = input2['u']
    input_rec2.T_g  = input2['T_g']
    input_rec2.S_n  = input2['S_n']

    m_pp         = input1['m_pp']
    percent_snow = input1['percent_snow']
    rho_snow     = input1['rho_snow']
    T_pp         = input1['T_pp']

    precip_now = 0
    if m_pp > 0:
        precip_now = 1

    # get the model state
    elevation       = output_rec['elevation']
    z_0             = output_rec['z_0']
    z_s             = output_rec['z_s']
    rho             = output_rec['rho']
    T_s_0           = output_rec['T_s_0']
    T_s_l           = output_rec['T_s_l']
    T_s             = output_rec['T_s']
    h2o_sat         = output_rec['h2o_sat']
    layer_count     = output_rec['layer_count']

    R_n_bar         = output_rec['R_n_bar']
    H_bar           = output_rec['H_bar']
    L_v_E_bar       = output_rec['L_v_E_bar']
    G_bar           = output_rec['G_bar']
    G_0_bar         = output_rec['G_0_bar']
    M_bar           = output_rec['M_bar']
    delta_Q_bar     = output_rec['delta_Q_bar']
    delta_Q_0_bar   = output_rec['delta_Q_0_bar']
    E_s_sum         = output_rec['E_s_sum']
    melt_sum        = output_rec['melt_sum']
    ro_pred_sum     = output_rec['ro_pred_sum']

#         print z_0

    # establish conditions for snowpack
    # the first step mimic's snobal which only calls init_snow once. This
    # might mean that the model states will need to be saved in memory
    # or there will be a slight descrepancy with isnobal. But with this,
    # there should be a descrepancy in isnobal as well
    if first_step:
        init_snow()

    # set air pressure from site elev
    P_a = HYSTAT(SEA_LEVEL, STD_AIRTMP, STD_LAPSE, (elevation / 1000.0),
        GRAVITY, MOL_AIR)

    # do_data_tstep.c
    # dt = do_data_tstep()
    dt = cython_do_data_tstep(input1, input2)
    if dt == 0:
        rt = False
        abort()

    output_rec['current_time'] = current_time
    output_rec['time_since_out'] = time_since_out

    output_rec['elevation'] = elevation
    output_rec['z_0'] = z_0
    output_rec['rho'] = rho
    output_rec['T_s_0'] = T_s_0
    output_rec['T_s_l'] = T_s_l
    output_rec['T_s'] = T_s
    output_rec['h2o_sat'] = h2o_sat
    output_rec['h2o_max'] = h2o_max
    output_rec['h2o'] = h2o
    output_rec['layer_count'] = layer_count
    output_rec['cc_s_0'] = cc_s_0
    output_rec['cc_s_l'] = cc_s_l
    output_rec['cc_s'] = cc_s
    output_rec['m_s_0'] = m_s_0
    output_rec['m_s_l'] = m_s_l
    output_rec['m_s'] = m_s
    output_rec['z_s_0'] = z_s_0
    output_rec['z_s_l'] = z_s_l
    output_rec['z_s'] = z_s

    output_rec['R_n_bar'] = R_n_bar
    output_rec['H_bar'] = H_bar
    output_rec['L_v_E_bar'] = L_v_E_bar
    output_rec['G_bar'] = G_bar
    output_rec['G_0_bar'] = G_0_bar
    output_rec['M_bar'] = M_bar
    output_rec['delta_Q_bar'] = delta_Q_bar
    output_rec['delta_Q_0_bar'] = delta_Q_0_bar
    output_rec['E_s_sum'] = E_s_sum
    output_rec['melt_sum'] = melt_sum
    output_rec['ro_pred_sum'] = ro_pred_sum

    return rt


def cython_do_data_tstep(input1, input2):
	# //	static PRECIP_REC *pp_info    = precip_info;
	# /* precip info for data timestep */
	# //	static TSTEP_REC  *data_tstep = tstep_info;
	# /* timestep info for data timestep */
	# static PRECIP_REC *pp_info
	# pp_info = precip_info

	# static TSTEP_REC *data_tstep
	# data_tstep = tstep_info

    cdef int level;			# loop index

    # Copy values from first input record into global variables.
    S_n  = input1['S_n']
    I_lw = input1['I_lw']
    T_a  = input1['T_a']
    e_a  = input1['e_a']
    u    = input1['u']
    T_g  = input1['T_g']
    
    # Compute deltas for the climate input parameters over
    # the data timestep.
    input_deltas[DATA_TSTEP].S_n  = input_rec2.S_n  - input_rec1.S_n
    input_deltas[DATA_TSTEP].I_lw = input_rec2.I_lw - input_rec1.I_lw
    input_deltas[DATA_TSTEP].T_a  = input_rec2.T_a  - input_rec1.T_a
    input_deltas[DATA_TSTEP].e_a  = input_rec2.e_a  - input_rec1.e_a
    input_deltas[DATA_TSTEP].u    = input_rec2.u    - input_rec1.u
    input_deltas[DATA_TSTEP].T_g  = input_rec2.T_g  - input_rec1.T_g

	# /*
	#  *  If there is precipitation, then compute the amount of rain &
	#  *  snow in it.
	#  */
	# if (precip_now) {
	# 	pp_info->m_pp   = m_pp;
	# 	pp_info->m_snow = percent_snow * m_pp;
	# 	pp_info->m_rain = m_pp - pp_info->m_snow;
	# 	if (pp_info->m_snow > 0.0) {
	# 		if (rho_snow > 0.0)
	# 			pp_info->z_snow = pp_info->m_snow / rho_snow;
	# 		else {
	# 			//				usrerr("rho_snow is <= 0.0 with %_snow > 0.0");
	# 			fprintf(stderr, "rho_snow is <= 0.0 with %_snow > 0.0");
	# 			return FALSE;
	# 		}
	# 	}
	# 	else
	# 		pp_info->z_snow = 0.0;

	# 	/*
	# 	 *  Mixed snow and rain
	# 	 */
	# 	if ((pp_info->m_snow > 0.0) && (pp_info->m_rain > 0.0)) {
	# 		T_snow = FREEZE;
	# 		h2o_sat_snow = 1.0;
	# 		T_rain = T_pp;
	# 	}

	# 	/*
	# 	 *  Snow only
	# 	 */
	# 	else if (pp_info->m_snow > 0.0) {
	# 		if (T_pp < FREEZE) {		/* Cold snow */
	# 			T_snow = T_pp;
	# 			h2o_sat_snow = 0.0;
	# 		}
	# 		else {				/* Warm snow */
	# 			T_snow = FREEZE;
	# 			h2o_sat_snow = 1.0;
	# 		}
	# 	}

	# 	/*
	# 	 *  Rain only
	# 	 */
	# 	else if (pp_info->m_rain > 0.0) {
	# 		T_rain = T_pp;
	# 	}
	# }

	# /*
	#  *  Clear the 'computed' flag at the other timestep levels.
	#  */
	# for (level = NORMAL_TSTEP; level <= SMALL_TSTEP; level++)
	# 	computed[level] = FALSE;

	# /*
	#  *  Divide the data timestep into normal run timesteps.
	#  */
	# return _divide_tstep(data_tstep);
