#ifndef	SOLAR_H
#define	SOLAR_H

/*
 *  Solar radiation and geometry.
 */

/* ------------------------------------------------------------------------ */

/*
 *  Maximum angle of declination (degrees).
 */
#define MAX_DECLINATION		23.6

/*
 *  solar constant from Willson's Eos presentation (in W / m^2)
 */
#define SOLAR_CONSTANT		1368.0

/* ------------------------------------------------------------------------ */
 
extern int		ephemeris(datetime_t *gmt_dt, double *r, double *declin,
				  double *omega);
extern int		sunpath(double lat, double lon, double declin,
			        double omega, double *cosZ, double *azm);
extern int		sunangle(double lat, double lon, datetime_t *dt,
				 double * mu, double * azmSun, double * radVec);
extern datetime_t *	sunrise(double lat, double lon, int year, int month,
				int day);
extern datetime_t *	sunset (double lat, double lon, int year, int month,
				int day);

#endif
