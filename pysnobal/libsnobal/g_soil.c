//#include "ipw.h"
#include "snow.h"

double
g_soil(
		double	rho,	/* snow layer's density (kg/m^3)	     */
		double	tsno,	/* snow layer's temperature (K)		     */
		double	tg,	/* soil temperature (K)			     */
		double	ds,	/* snow layer's thickness (m)		     */
		double	dg,	/* dpeth of soil temperature measurement (m) */
		double	pa)	/* air pressure (Pa)			     */
{
	double	k_g;
	double	kcs;
	double	k_s;
	double	g;

	/*	check tsno	*/
	if (tsno > FREEZE) {
//		warn("g_soil: tsno = %8.2f; set to %8.2f\n", tsno, FREEZE);
		fprintf(stdout, "g_soil: tsno = %8.2f; set to %8.2f\n", tsno, FREEZE);
		tsno = FREEZE;
	}

	/*	set effective soil conductivity	*/
	/***	changed to KT_MOISTSAND by D. Marks, NWRC, 09/30/2003	***/
	/***	based on heat flux data from RMSP			***/
	/***	note: Kt should be passed as an argument		***/
	/***	k_g = efcon(KT_WETSAND, tg, pa);			***/
	k_g = efcon(KT_MOISTSAND, tg, pa);

	/*	calculate G	*/
	/*	set snow conductivity	*/
	kcs = KTS(rho);
	k_s = efcon(kcs, tsno, pa);

	g = ssxfr(k_s, k_g, tsno, tg, ds, dg);

	return (g);
}
