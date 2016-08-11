#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <omp.h>
#include "envphys.h"

#define AH		1.0	/* ratio sensible/momentum phi func	*/
#define AV		1.0	/* ratio latent/momentum phi func	*/
#define ITMAX		50	/* max # iterations allowed		*/
#define PAESCHKE	7.35	/* Paeschke's const (eq. 5.3)		*/
#define THRESH		1.e-5	/* convergence threshold		*/

#define SM		0
#define SH		1
#define SV		2
#define BETA_S		5.2
#define BETA_U		16

/* ----------------------------------------------------------------------- */

/*
 * psi-functions
 *	code =	SM	momentum
 *		SH	sensible heat flux
 *		SV	latent heat flux
 */

double
psi(
		double	zeta,		/* z/lo				*/
		int	code)		/* which psi function? (see above) */
{
	double	x;		/* height function variable	*/
	double	result;

	if (zeta > 0) {		/* stable */
		if (zeta > 1)
			zeta = 1;
		result = -BETA_S * zeta;
	}

	else if (zeta < 0) {	/* unstable */

		x = sqrt(sqrt(1 - BETA_U * zeta));

		switch (code) {
		case SM:
			result = 2 * log((1+x)/2) + log((1+x*x)/2) -
			2 * atan(x) + M_PI_2;
			break;

		case SH:
		case SV:
			result = 2 * log((1+x*x)/2);
			break;

		default: /* shouldn't reach */
			perror("psi-function code not of these: SM, SH, SV");
			result = -1;
		}
	}

	else {			/* neutral */
		result = 0;
	}

	return (result);
}

/* ----------------------------------------------------------------------- */

int
hle1(
		double	press,	/* air pressure (Pa)			*/
		double	ta,	/* air temperature (K) at height za	*/
		double	ts,	/* surface temperature (K)		*/
		double	za,	/* height of air temp measurement (m)	*/
		double	ea,	/* vapor pressure (Pa) at height zq	*/
		double	es,	/* vapor pressure (Pa) at surface	*/
		double	zq,	/* height of spec hum measurement (m)	*/
		double	u,	/* wind speed (m/s) at height zu	*/
		double	zu,	/* height of wind speed measurement (m)	*/
		double	z0,	/* roughness length (m)			*/
		int error_check, /* flag if the inputs should be checked for errors */

		/* output variables */

		double *h,	/* sens heat flux (+ to surf) (W/m^2)	*/
		double *le,	/* latent heat flux (+ to surf) (W/m^2)	*/
		double *e)	/* mass flux (+ to surf) (kg/m^2/s)	*/
{
	double	ah = AH;
	double	av = AV;
	double	cp = CP_AIR;
	double	d0;	/* displacement height (eq. 5.3)	*/
	double	dens;	/* air density				*/
	double	diff;	/* difference between guesses		*/
	double	factor;
	double	g = GRAVITY;
	double	k = VON_KARMAN;
	double	last;	/* last guess at lo			*/
	double	lo;	/* Obukhov stability length (eq. 4.25)	*/
	double	ltsh;	/* log ((za-d0)/z0)			*/
	double	ltsm;	/* log ((zu-d0)/z0)			*/
	double	ltsv;	/* log ((zq-d0)/z0)			*/
	double	qa;	/* specific humidity at height zq	*/
	double	qs;	/* specific humidity at surface		*/
	double	ustar;	/* friction velocity (eq. 4.34')	*/
	double	xlh;	/* latent heat of vap/subl		*/
	int	ier;	/* return error code			*/
	int	iter;	/* iteration counter			*/
	double e_ts, e_ta;

	/*
	 * check for bad input
	 */

	if (error_check == 1) {

		/* heights must be positive */
		if (z0 <= 0 || zq <= z0 || zu <= z0 || za <= z0) {
			perror (printf("height not positive; z0=%f\tzq=%f\tzu=%\tza=%f",
					z0, zq, zu, za));
			ier = -2;
			return (ier);
		}

		/* temperatures are Kelvin */
		if (ta <= 0 || ts <= 0) {
			perror (printf("temps not K; ta=%f\tts=%f", ta, ts));
			ier = -2;
			return (ier);
		}

		/* pressures must be positive */
		if (ea <= 0 || es <= 0 || press <= 0 || ea >= press || es >= press) {
			perror (printf("press < 0; ea=%f\tes=%f\tpress=%f", ea, es, press));
			ier = -2;
			return (ier);
		}

		/* vapor pressures can't exceed saturation */
		/* if way off stop */
		e_ts = sati(ts);
		e_ta = satw(ta);
		if ((es - 25.0) > e_ts || (ea - 25.0) > e_ta) {
			perror (printf("vp > sat; es=%f\tessat=%f\tea=%f\teasat=%f",
					es, e_ts, ea, sati(ta)));
			ier = -2;
			return (ier);
		}

		/* else fix them up */
		if (es > e_ts) {
			es = e_ts;
		}
		if (ea > e_ta) {
			ea = e_ta;
		}
	}

	/*
	 * displacement plane height, eq. 5.3 & 5.4
	 */

	d0 = 2 * PAESCHKE * z0 / 3;

	/*
	 * constant log expressions
	 */

	ltsm = log((zu - d0) / z0);
	ltsh = log((za - d0) / z0);
	ltsv = log((zq - d0) / z0);

	/*
	 * convert vapor pressures to specific humidities
	 */
	qa = SPEC_HUM(ea, press);
	qs = SPEC_HUM(es, press);

	/*
	 * convert temperature to potential temperature
	 */

	ta += DALR * za;

	/*
	 * air density at press, virtual temp of geometric mean
	 * of air and surface
	 */

	dens = GAS_DEN(press, MOL_AIR,
			VIR_TEMP(sqrt(ta*ts), sqrt(ea*es), press));

	/*
	 * starting value, assume neutral stability, so psi-functions
	 * are all zero
	 */

	ustar = k * u / ltsm;
	factor = k * ustar * dens;
	*e = (qa - qs) * factor * av / ltsv;
	*h = (ta - ts) * factor * cp * ah / ltsh;

	/*
	 * if not neutral stability, iterate on Obukhov stability
	 * length to find solution
	 */

	iter = 0;
	if (ta != ts) {

		lo = HUGE_VAL;

		do {
			last = lo;

			/*
			 * Eq 4.25, but no minus sign as we define
			 * positive H as toward surface
			 */

			/*
			 * There was an error in the old version of this
			 * line that omitted the cubic power of ustar.
			 * Now, this error has been fixed.
			 */

			lo = ustar * ustar * ustar * dens 
					/ (k * g * (*h/(ta*cp) + 0.61 * *e));

			/*
			 * friction velocity, eq. 4.34'
			 */

			ustar = k * u / (ltsm - psi(zu/lo, SM));

			/*
			 * evaporative flux, eq. 4.33'
			 */

			factor = k * ustar * dens;
			*e = (qa - qs) * factor * av /
					(ltsv - psi(zq/lo, SV));

			/*
			 * sensible heat flus, eq. 4.35'
			 * with sign reversed
			 */

			*h = (ta - ts) * factor * ah * cp /
					(ltsh - psi(za/lo, SH));

			diff = last - lo;

		} while (fabs(diff) > THRESH &&
				fabs(diff/lo) > THRESH &&
				++iter < ITMAX);
	}

	ier = (iter >= ITMAX)? -1 : 0;

	xlh = LH_VAP(ts);
	if (ts <= FREEZE)
		xlh += LH_FUS(ts);

	/*
	 * latent heat flux (- away from surf)
	 */
	*le = xlh * *e;

	return (ier);
}


/* ----------------------------------------------------------------------- */

int
hle1_grid(
		int ngrid,		/* number of points on grid */
		double	*press,	/* air pressure (Pa)			*/
		double	*ta,	/* air temperature (K) at height za	*/
		double	*ts,	/* surface temperature (K)		*/
		double	*za,	/* height of air temp measurement (m)	*/
		double	*ea,	/* vapor pressure (Pa) at height zq	*/
		double	*es,	/* vapor pressure (Pa) at surface	*/
		double	*zq,	/* height of spec hum measurement (m)	*/
		double	*u,	/* wind speed (m/s) at height zu	*/
		double	*zu,	/* height of wind speed measurement (m)	*/
		double	*z0,	/* roughness length (m)			*/
		int error_check, /* flag if the inputs should be checked for errors */
		int nthreads, 	/* number of threads to use */

		/* output variables */

		double *h,	/* sens heat flux (+ to surf) (W/m^2)	*/
		double *le,	/* latent heat flux (+ to surf) (W/m^2)	*/
		double *e)	/* mass flux (+ to surf) (kg/m^2/s)	*/
{


	int i, ier;
	int brokeit = 0;

	// setup the parallel execution
	omp_set_dynamic(100);     // Explicitly disable dynamic teams
	omp_set_num_threads(nthreads); // Use N threads for all consecutive parallel regions

#pragma omp parallel shared(press, ta, ts, za, ea, es, zq, u, zu, z0, error_check, h, le, e, brokeit) private(i, ier)
	{
#pragma omp for
		for (i =0; i < ngrid; i++) {

			ier = hle1(press[i], ta[i], ts[i], za[i], ea[i], es[i], zq[i], u[i], zu[i], z0[i], error_check, &h[i], &le[i], &e[i]);

			if (ier == -1) {
				perror(printf("hle1 did not converge at point %i", i));
				brokeit = i;
//				return i;
			}

		}
	}

	ier = brokeit;

	return ier;

}

/* ----------------------------------------------------------------------- */

double sati_grid(
		int ngrid, 		/* number of points */
		double *tk,		/* temperature (K) */
		double *es		/* saturation vapor pressure output (Pa) */
)
{

	int i;
	double x;

	for (i =0; i < ngrid; i++) {

		x = sati(tk[i]);

		if (x == -1) {
			perror(printf("hle1 did not converge at point %i", i));
			return i;
		}
		es[i] = x;

	}

	return x;

}



double
sati(
		double  tk)		/* air temperature (K)	*/
{
	double  l10;
	double  x;

	if (tk <= 0.) {
		printf("tk=%f\n",tk);
		assert(tk > 0.);
	}

	if (tk > FREEZE) {
		x = satw(tk);
		return(x);
	}

	errno = 0;
	l10 = log(1.e1);

	x = pow(1.e1,-9.09718*((FREEZE/tk)-1.) - 3.56654*log(FREEZE/tk)/l10 +
			8.76793e-1*(1.-(tk/FREEZE)) + log(6.1071)/l10);

	if (errno) {
		perror("sati: bad return from log or pow");
		return (-1);
	}

	return(x*1.e2);
}


double
satw(
		double  tk)		/* air temperature (K)		*/
{
	double  x;
	double  l10;

	if (tk <= 0.) {
		assert(tk > 0.);
	}

	errno = 0;
	l10 = log(1.e1);

	x = -7.90298*(BOIL/tk-1.) + 5.02808*log(BOIL/tk)/l10 -
			1.3816e-7*(pow(1.e1,1.1344e1*(1.-tk/BOIL))-1.) +
			8.1328e-3*(pow(1.e1,-3.49149*(BOIL/tk-1.))-1.) +
			log(SEA_LEVEL)/l10;

	x = pow(1.e1,x);

	if (errno) {
		perror("satw: bad return from log or pow");
		return (-1);
	}

	return(x);
}


/* ----------------------------------------------------------------------- */

void efcon_grid (
		int ngrid, 		/* number of points */
		double *k,		/* layer thermal conductivity (J/(m K sec)) */
		double *t,		/* layer temperature (K)		    */
		double *p,		/* air pressure (Pa)  			    */
		double *e,		/* layer vapor pressure (Pa)		*/
		double *etc 	/* output, effecitve layer diffusion */
)
{
	/*	calculate effective layer diffusion
		(see Anderson, 1976, pg. 32)		*/

	int i;
	double	de;
	double	lh;
	double	q;

	for (i = 0; i < ngrid; i++) {
		de = DIFFUS(p[i], t[i]);

		/*	set latent heat from layer temp.	*/
		if(t[i] > FREEZE)
			lh = LH_VAP(t[i]);
		else if(t[i] == FREEZE)
			lh = (LH_VAP(t[i]) + LH_SUB(t[i])) / 2.0;
		else
			lh = LH_SUB(t[i]);

		/*	set mixing ratio from layer temp.	*/
		q = MIX_RATIO(e[i], p[i]);

		/*	calculate effective layer conductivity	*/
		etc[i] = k[i] + (lh * de * q);
	}


}
