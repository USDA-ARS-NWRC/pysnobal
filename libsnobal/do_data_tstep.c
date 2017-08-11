/*
 ** NAME
   do_data_tstep -- run model for 1 data timestep between 2 input records
 **
 ** SYNOPSIS
 **	#include "snobal.h"
 **
 **	int
 **	do_data_tstep(void)
 **
 ** DESCRIPTION
 **	This routine performs the model's calculations for 1 data timestep
 **	between 2 input-data records which are in 'input_rec1' and
 **	'input_rec2'.
 **
 **	If there's precipitation during the data timestep, the flag
 **	'precip_now' used be TRUE.  Furthermore, the routine requires
 **	that the following precipitation variables have been initialized:
 **
 **		m_pp
 **		percent_snow
 **		rho_snow
 **		T_pp
 **
 **	This routine divides the data timestep into the appropriate number
 **	of normal run timesteps.  The input values for each normal timestep
 **	are computed from the two input records by linear interpolation.
 **
 **	If output is desired for any of the run timesteps (normal, medium,
 **	or small), the appropriate output flags must be set in the proper
 **	timestep's record (i.e., the array 'tstep_info').  If any output
 **	flag is set, the routine requires that the global variable 'out_func'
 **	point to appropriate output function.
 **
 **	This routine may return in the middle of a data timestep if:
 **
 **		a)  the output function pointed to by 'out_func' is called, and
 **		b)  the flag 'run_no_snow' is FALSE, and
 **		c)  there is no snow remaining on the ground at the end of
 **		    timestep
 **
 **	In this happens, the flag 'stop_no_snow' is set to TRUE.
 **
 ** RETURN VALUE
 **
 **	TRUE	The model's calculations were completed.
 **
 **	FALSE	An error occured, and a message explaining the error has
 **		been stored with the 'usrerr' routine.
 **
 ** GLOBAL VARIABLES READ
 **	e_a
 **	I_lw
 **	in_rec
 **	layer_count
 **	m_pp_data
 **	m_rain_data
 **	m_snow_data
 **	more_pr_recs
 **	precip_data
 **	ro
 **	ro_data
 **	run_no_snow
 **	S_n
 **	T_a
 **	T_g
 **	tstep_info
 **	u
 **	z_snow_data
 **
 ** GLOBAL VARIABLES MODIFIED
 **	precip_now
 **	stop_no_snow
 */

#include <omp.h>
//#include	"ipw.h"
#include        "envphys.h"
#include        "_snobal.h"

int
do_data_tstep(void)
{
	//	static PRECIP_REC *pp_info    = precip_info;
	/* precip info for data timestep */
	//	static TSTEP_REC  *data_tstep = tstep_info;
	/* timestep info for data timestep */
	static PRECIP_REC *pp_info;
#pragma omp threadprivate(pp_info)
	pp_info = precip_info;

	static TSTEP_REC *data_tstep;
#pragma omp threadprivate(data_tstep)
	data_tstep = tstep_info;

//	printf("%i -- %i -- %f -- %f\n", tstep_info[0].level, tstep_info[0].time_step, tstep_info[0].intervals, tstep_info[0].threshold);
//	printf("%i -- %i -- %f -- %f\n", tstep_info[1].level, tstep_info[1].time_step, tstep_info[1].intervals, tstep_info[1].threshold);
//	printf("%i -- %i -- %f -- %f\n", tstep_info[2].level, tstep_info[2].time_step, tstep_info[2].intervals, tstep_info[2].threshold);
//	printf("%i -- %i -- %f -- %f\n", tstep_info[3].level, tstep_info[3].time_step, tstep_info[3].intervals, tstep_info[3].threshold);
//
//	printf("%i -- %i -- %f -- %f\n", data_tstep[0].level, data_tstep[0].time_step, data_tstep[0].intervals, data_tstep[0].threshold);
//	printf("%i -- %i -- %f -- %f\n", data_tstep[1].level, data_tstep[1].time_step, data_tstep[1].intervals, data_tstep[1].threshold);
//	printf("%i -- %i -- %f -- %f\n", data_tstep[2].level, data_tstep[2].time_step, data_tstep[2].intervals, data_tstep[2].threshold);
//	printf("%i -- %i -- %f -- %f\n", data_tstep[3].level, data_tstep[3].time_step, data_tstep[3].intervals, data_tstep[3].threshold);

	int	level;			/* loop index */

	/*
	 *  Copy values from first input record into global variables.
	 */
	S_n  = input_rec1.S_n;
	I_lw = input_rec1.I_lw;
	T_a  = input_rec1.T_a;
	e_a  = input_rec1.e_a;
	u    = input_rec1.u  ;
	T_g  = input_rec1.T_g;
	if (ro_data)
		ro = input_rec1.ro;

	//	printf("%f - %f - %f - %f - %f - %f\n", S_n, I_lw, T_a, e_a, u, T_g);
	//	printf("%f\n", max_h2o_vol);
	//	printf("%i %i %i %i\n", DATA_TSTEP, NORMAL_TSTEP, MEDIUM_TSTEP, SMALL_TSTEP);

	/*
	 *  Compute deltas for the climate input parameters over
	 *  the data timestep.
	 */
	input_deltas[DATA_TSTEP].S_n  = input_rec2.S_n  - input_rec1.S_n;
	input_deltas[DATA_TSTEP].I_lw = input_rec2.I_lw - input_rec1.I_lw;
	input_deltas[DATA_TSTEP].T_a  = input_rec2.T_a  - input_rec1.T_a;
	input_deltas[DATA_TSTEP].e_a  = input_rec2.e_a  - input_rec1.e_a;
	input_deltas[DATA_TSTEP].u    = input_rec2.u    - input_rec1.u;
	input_deltas[DATA_TSTEP].T_g  = input_rec2.T_g  - input_rec1.T_g;
	if (ro_data)
		input_deltas[DATA_TSTEP].ro = input_rec2.ro - input_rec1.ro;

	/*
	 *  If there is precipitation, then compute the amount of rain &
	 *  snow in it.
	 */
	if (precip_now) {
		pp_info->m_pp   = m_pp;
		pp_info->m_snow = percent_snow * m_pp;
		pp_info->m_rain = m_pp - pp_info->m_snow;
		if (pp_info->m_snow > 0.0) {
			if (rho_snow > 0.0)
				pp_info->z_snow = pp_info->m_snow / rho_snow;
			else {
				//				usrerr("rho_snow is <= 0.0 with %_snow > 0.0");
				fprintf(stderr, "rho_snow is <= 0.0 with %_snow > 0.0");
				return FALSE;
			}
		}
		else
			pp_info->z_snow = 0.0;

		/*
		 *  Mixed snow and rain
		 */
		if ((pp_info->m_snow > 0.0) && (pp_info->m_rain > 0.0)) {
			T_snow = FREEZE;
			h2o_sat_snow = 1.0;
			T_rain = T_pp;
		}

		/*
		 *  Snow only
		 */
		else if (pp_info->m_snow > 0.0) {
			if (T_pp < FREEZE) {		/* Cold snow */
				T_snow = T_pp;
				h2o_sat_snow = 0.0;
			}
			else {				/* Warm snow */
				T_snow = FREEZE;
				h2o_sat_snow = 1.0;
			}
		}

		/*
		 *  Rain only
		 */
		else if (pp_info->m_rain > 0.0) {
			T_rain = T_pp;
		}
	}

	/*
	 *  Clear the 'computed' flag at the other timestep levels.
	 */
	for (level = NORMAL_TSTEP; level <= SMALL_TSTEP; level++)
		computed[level] = FALSE;

	/*
	 *  Divide the data timestep into normal run timesteps.
	 */
	return _divide_tstep(data_tstep);
}
