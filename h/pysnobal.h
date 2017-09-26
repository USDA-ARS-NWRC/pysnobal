/*
 ** NAME
 **      pgm.h
 **
 ** DESCRIPTION
 **      The include file for 'isnobal'.
 */

#ifndef _ISNOBAL_H_
#define _ISNOBAL_H_

#define DEFAULT_Z_U	5.0	/* default wind speed measurement height */
#define DEFAULT_Z_T	5.0	/* default air temp and vapor press hght */

#define IBANDS		6	/* # bands in input image	 	 */
#define EMBANDS	       10	/* # bands in energy/mass output image	 */
#define SBANDS	        9	/* # bands in snow output image		 */
#define PBANDS		4	/* # bands in precip image		 */
#define ICBANDS		7	/* # bands in initial conditions image   */
#define ICBANDS_RESTART	8	/* # bands in init cond image (restart)  */
#define TBANDS	       17	/* # bands in temporary results file	 */

#define NO_DATA   -999999	/* output value for masked pnt (no data) */

typedef struct {
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
} OUTPUT_REC;
//typedef OUTPUT_REC *out_p;
//extern OUTPUT_REC output_rec[100];	/* output data structure */

typedef struct {
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
	double* h2o_vol;
	double* h2o;
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
} OUTPUT_REC_ARR;

typedef struct {
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
} INPUT_REC_ARR;


typedef struct {
	double z_u;
	double z_T;
	double z_g;
	int relative_heights;
	double max_h2o_vol;
	double max_z_s_0;
} PARAMS;

/* ------------------------------------------------------------------------- */

/*
 *  Routines that are part of isnobal program.
 */

//extern int call_snobal(int N, int nthreads, int first_step, TSTEP_REC tstep_info[4], OUTPUT_REC** output_rec, INPUT_REC_ARR* input1, INPUT_REC_ARR* input2, PARAMS params, OUTPUT_REC_ARR* output1);
extern int call_snobal(int N, int nthreads, int first_step, TSTEP_REC tstep_info[4], INPUT_REC_ARR* input1, INPUT_REC_ARR* input2, PARAMS params, OUTPUT_REC_ARR* output1);

//extern	void	assign_buffers (int masked, int n, int output, OUTPUT_REC **output_rec);
//extern	void	buffers        (void);
//extern	void	check_range    (int index, double value, double min, double max,
//		char * descrip, bool_t print_line_samp);
//extern	void	check_units    (LQH_T **lq_headers, UNITS_T *units, int nbands,
//		int fd);
//extern	void	copy_image     (char *tempfile, int nbands, fpixel_t * buf,
//		int fdo);
//extern	void 	e_m_image      (int step, OUTPUT_REC **output_rec, int nbits);
//extern	bool_t	extract_data   (bool_t first_step, int n, bool_t sun_up[], OUTPUT_REC **output_rec);
//extern	void	headers        (void);
//extern	void	isnobal		   (int out_step, int nthreads, int dynamic_teams, int got_opt_F, int verbose, int nbits);
///*extern	void	isnobal        (int out_step);*/
//extern	void	newlqh         (int fdo, int nbands, fpixel_t *mins,
//		fpixel_t *maxs, char **units);
//extern	int	open_input     (char *prefix, int index, bool_t *sun_up);
//extern	int	output_image   (char * filename, int nbands, char ** units,
//		char ** annots, fpixel_t * mins,
//		fpixel_t * maxs, int nbits);
//extern	bool_t	precip_event   (float curr_time, char *pre_img);
//extern	void	precip_hdrs    (char *filename);
//extern	void	read_data      (int first_step);
//extern	void 	snow_image     (int step, OUTPUT_REC **output_rec, int nbits);
//extern	void	temp_filename  (char *prefix, char *filename);
//extern	void	warn_range    (int index, double value, double min, double max,
//		char * descrip, bool_t print_line_samp);
//extern	void	write_data     (int output, int last_step);

/* ------------------------------------------------------------------------- */

/*
 *  Global variables internal to isnobal program.
 */

extern	int	 	units_warn;	/* check units in input images?     */
extern	char	       *compress_cmd;	/* shell command to compress images */

/* timesteps and indices */

extern	int		start_step;	/* index of first timestep	     */
extern	int		nstep;		/* #  of data timesteps		     */
extern	int		nDigits;	/* #  of digits in suffixes of images*/
extern	bool_t		restart;	/* restart flag			     */

/* model variables */


extern double elevation;
//#pragma omp threadprivate(elevation)

#endif /* _ISNOBAL_H_ */
