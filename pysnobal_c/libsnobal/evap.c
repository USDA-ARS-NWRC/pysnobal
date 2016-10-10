//#include "ipw.h"
#include "envphys.h"

double
evap(
	double	le,	/* latent heat transfer (W/m**2) */
	double	ts)	/* surface temperature (K)       */
{
	double	rate;
	double	lh;

	if(ts > FREEZE)
		lh = LH_VAP(ts);
	else if(ts == FREEZE)
		lh = (LH_VAP(ts) + LH_SUB(ts)) / 2.0;
	else
		lh = LH_SUB(ts);

	rate = le / lh;

	return(rate);
}
