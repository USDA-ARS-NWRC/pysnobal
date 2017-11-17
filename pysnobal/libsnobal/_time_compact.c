/*
** NAME
**      _time_compact -- snowcover gravety, depth & temperature compaction
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_time_compact(void)
**
** DESCRIPTION
**	This routine replaces the original simple gravety compaction routine
**	which "aged" the snowcover by accounting for the compaction
**	or densification by gravity as time passes.  The original routine
**	relied on an approximation based on equations in Anderson (1976), which
**	increased the snowcover's density using the following "half-saturation"
**	function:
**
**		rho(time) = A / (1 + B/time)
**
**	Historically, precise snow density was not a concern, as long as mass
**	and SWE were correct.  With the development of the ASO program providing
**	time-series lidar snow depth images, and the 2017 snow season in the
*8	southern Sierra Nevada, with individual storms approaching 500mm of
**	deposition in a few days, and upper elevation snow depths of greater
**	than 10m, it became clear that a more robust density model was required.
**
**	Snow Density:  Snow density is initially a function of the temperature
**	of the ice particles (snow flakes) as they fall to the ground during a
**	storm.  Under very cold conditions (-10 to -15 C), we can get snow as
**	light as 50 kg m-3, which is really light powder and great skiing.
**	As the ice particle temperature approaches 0C, new snow densities can
**	be as high as 200 kg m-3, which is not light powder, but still pretty
**	good skiing, unless there is a lot of it. There are several (4)
**	processes that can increase snow density during and after the storm.
**	Note that the largest and most rapid changes occur during or just
**	following the storm.  Compaction (an increase in density without a
**	change in mass) can be caused by:
**
**	1) Destructive mechanical metamorphism (compaction due to wind -
**	   affects mainly low density near-surface or new snow)
**	2) Destructive temperature metamorphism
**	3) Pressure - compaction due to snow load or overburden (affects both
**	   new snow deposition, snow on the ground, including drifting)
**	4) Addition of liquid water due to melt or rain
**
**	Of these compaction processes, iSnobal accounts three - 2,3 & 4.
**	This routine addresses 2, temperature metamorphism, and 3., overburden
**	metamorphism. The 4th, addition of liquid water is accounted for in the
**	routine _h2o_compact.c. We are using equations found in Anderson (1976)
**	and Oleson, et al. (2013), who got his equations from Anderson in the
**	first place.  (Oleson made it easier to figure out the units...)
**
**	Though many dedicated individuals have worked on the density issue over
**	over the last 60+ years, we have, in general, a group working in Saporro
**	Japan to thank for most of what we know about snow density. Anderson
**	got the data and basic fitting equations from careful field and cold
**	room measurements made by Yosida (1963), Mellor (1964) and Kojima (1967).
**	The equations we are using are based on those that they derived from data
**	that they carefully collected over a number of years.  It is noteworthy
**	that while rapidly changing computers have made the kind of spatial
**	modeling we are attempting possible, snow physics remains unchanged – and
**	the field and labortory efforts and data fitting equations from more than
**	half a century ago represent our best understanding of those physics.
**
**	Tz = 0.0 (freezing temperature, C or K)
**	Ts = snow or precipitation temperature (C or K)
**	rho = intital snow density (kg/(m^3))
**	SWE = snow mass (mm/(m^2))
**	K = temperature metamorphism coef.
**	rho_n = new snow density (kg/(m^3))
**	zs = new snow depth (mm)
**
**	Proportional Destructive Temperature Metamorphism (PTM):
**
**	if (rho < 100)
**		K = 1.0
**	else
**		K = exp(-0.046 * (rho - 100))
**
**	PTM = 0.01 * K * exp(-0.04 * (Tz - Ts))
**
**	Proportional Overburden Compaction (POC):
**
**	POC = (0.026 * exp(-0.08 * (Tz - Ts)) * SWE * exp(-21.0 * rho))
**
**	New snow density and depth
**
**	rho_n = rho + ((PTM + POC) * rho)
**	zs_n = SWE / rho_n
**
** GLOBAL VARIABLES READ
**	time_step, T_s, m_s, rho
**
** GLOBAL VARIABLES MODIFIED
**	rho
**
*/

//#include        "ipw.h"
#include        "_snobal.h"
#include	"envphys.h"

#define PI	3.14159265

#define	RMX	550
	/*
	 *  Maximum density due to compaction (kg/m^3).
	 */

#define	R	48
#define R1	23.5
#define R2	24.5
	/*
	 *  days over which compaction occurs
	 */
#define SWE_MAX	2000.0
	/*
	 *  max swe to consider - after this swe value, compaction is maximized
	 */

#define water	1000.0
	/*
	 *  Density of water (kg/m^3).
	 */

#define hour	3600
	/*
	 *  seconds in an hour
	 */
void
_time_compact(void)
{
	double	c11;	/* temperature metamorphism coefficient (Anderson, 1976) */
	double	Tz;	/* Freezing temperature (K) */
	double	d_rho_m;
	double	d_rho_c;
	double	rate;

	/*
	 *  If the snow is already at or above the maximum density due to
	 *  compaction, then just leave.
	 */
	if ((!snowcover) || (rho >= RMX))
		return;

	Tz = FREEZE;

	/*
	 *  Calculate rate which compaction will be applied per time step.
	 *  Rate will be adjusted as time step varies.
	 */
	if (m_s >= SWE_MAX)
		rate = 1.0;
	else {
		rate = R1 * cos((PI * m_s) / SWE_MAX) + R2;
		rate = rate / (time_step / hour);
	}

	/** Proportional Destructive Temperature Metamorphism (d_rho_m) **/

	if (rho < 100)
		c11 = 1.0;
	else
		c11 = exp(-0.046 * (rho - 100));

	d_rho_m = 0.01 * c11 * exp(-0.04 * (Tz - T_s));
	d_rho_m /= rate;

	/** Proportional Overburden Compaction (d_rho_c) **/

	d_rho_c = (0.026 * exp(-0.08 * (Tz - T_s)) * m_s * exp(-21.0 * (rho / water)));
	d_rho_c /= rate;

	/**	Compute New snow density	**/

	rho = rho + ((d_rho_m + d_rho_c) * rho);

        /*
	 *  Adjust the snowcover for this new density.
	 */
	_new_density();
}
