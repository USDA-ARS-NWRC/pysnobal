/*
 * call_snobal.c takes input from a Python function and calls all the necessary C functions
 */

#include <stdio.h>

#include <math.h>
#include <string.h>
#include <omp.h>
#include "snobal.h"
#include "pysnobal.h"

int call_snobal (
		int N,
		TSTEP_REC* tstep_info,
		OUTPUT_REC** output_rec,
		INPUT_REC_ARR* input1
)
{
	int n;
	int nthreads = 1;

	double current_time, time_since_out;

	// Some debuging stuff
//	printf("%f - %f - %f - %f - %f - %f\n", input1->S_n[0], input1->I_lw[0], input1->T_a[0], input1->e_a[0], input1->u[0], input1->T_g[0]);
	//	printf ("%i -- %f\n", N, output_rec[0]->elevation);

	/* set threads */
	if (nthreads != 1) {
		omp_set_num_threads(nthreads); 	// Use N threads for all consecutive parallel regions
	}



	//#pragma omp parallel shared(output_rec, ibuf1, ibuf2, icbuf, pbuf, mbuf, sbuf, embuf, fdp, fdm, restart, output, first_step)\
	//		private(n) \
	//		copyin(tstep_info, z_u, z_T, z_g, relative_hts, max_z_s_0, max_h2o_vol, out_func)
	//			{
	//#pragma omp for schedule(dynamic, dynamic_teams)
	for (n = 0; n < N; n++) {

		/* initialize some global variables for
			   'snobal' library for each pass since
			   the routine 'do_data_tstep' modifies them */

		current_time = output_rec[n]->current_time;
		time_since_out = output_rec[n]->time_since_out;

		// get the input records
//		input_rec1.I_lw = input1['I_lw'][i,j];
//		input_rec1.T_a  = input1['T_a'][i,j];
//		input_rec1.e_a  = input1['e_a'][i,j];
//		input_rec1.u    = input1['u'][i,j];
//		input_rec1.T_g  = input1['T_g'][i,j];
//		input_rec1.S_n  = input1['S_n'][i,j];
//
//		input_rec2.I_lw = input2['I_lw'][i,j];
//		input_rec2.T_a  = input2['T_a'][i,j];
//		input_rec2.e_a  = input2['e_a'][i,j];
//		input_rec2.u    = input2['u'][i,j];
//		input_rec2.T_g  = input2['T_g'][i,j];
//		input_rec2.S_n  = input2['S_n'][i,j];
//
//		m_pp         = input1['m_pp'][i,j];
//		percent_snow = input1['percent_snow'][i,j];
//		rho_snow     = input1['rho_snow'][i,j];
//		T_pp         = input1['T_pp'][i,j];



//		precip_now = (fdp != ERROR);

		/* extract data from I/O buffers */

//		if (extract_data(first_step, n, sun_up, output_rec)) {
//
//			/* run model on data for this pixel */
//
//			if (! do_data_tstep())
//				error("During step %d, at pixel %i", step, n);
//
//			/* assign data to output buffers */
//
//			assign_buffers(FALSE, n, output, output_rec);
//
//		} else { /* masked point */
//			assign_buffers(TRUE, n, output, output_rec);
//		}

	}  /* for loop on grid */
	//			}







	return 1;

}
