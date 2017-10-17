/*
** NAME
**      _vars.c - private global variables for snow-balance library
**
** DESCRIPTION
**      Defines and allocates space for variables declared in "_snobal.h"
**
** HISTORY
**      July, 1985:  	written by D. Marks, CSL, UCSB;
**	January, 1990:  added variables necessary for new input files
**			by Kelly Longley, OSU;
**	May 1995	Converted to IPW by J.Domingo, OSU
**	Feb 1996	Changed variable names to match notation in
**			snowmelt papers (e.g. Marks and Dozier, 1992).
**			J. Domingo, OSU
**	Jun 1996	Separated into public and private variables for
**			snobal library.  J. Domingo, OSU
**	Apr 2008	Initialized one variable (the last one) so that all
**			variables are visible to linker on OS X.
**			J. Domingo, Green Code
*/

//#include        "ipw.h"
#include        "_snobal.h"

/* --------------------------------------------------------------------- */

/*
 *  Global variables that are private (internal) to the snobal library.
 */

	INPUT_REC input_deltas[4];	/* deltas for climate-input parameters
					   over each timestep */

	PRECIP_REC precip_info[4];	/* array of precip info adjusted for
					   each timestep */

	int  computed[4];		/* array of flags for each timestep;
					   TRUE if computed values for input
					   deltas and precip arrays */

//	double  h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */
	int	isothermal;	/* melting? */
        int     snowcover = 0;	/* snow on gnd at start of current timestep?
                                     (see below about initialization) */

/*
 * One variable has been selected in this file and the "vars.c" file for
 * initialization to solve a linker problem on Mac OS X.  The "ar"
 * utility does not put these global variables in the library's symbol
 * table (index).  Therefore, the linker cannot find them when linking
 * programs like snobal and isnobal to this library.  This situation is
 * described in this message on GNU libtool mailing list:
 *
 *   "Mac OS X static linking (ranlib) fix"
 *   http://osdir.com/ml/gnu.libtool.general/2002-07/msg00045.html
 *
 * It is also described on this page about compilation problems at the
 * Fink project's web site (see the problem titled "I'm getting undefined
 * symbols but they are defined dammit"):
 *
 *   http://wiki.finkproject.org/index.php/Metapkg:compilation_problems
 *
 * One of the recommended solutions -- ranlib's "-c" option -- does solve
 * the linker problem.  But ranlib is not otherwise needed to build IPW.
 * The documentation of another solution -- gcc's "-fno-common" option --
 * lead to the current (hopefully long-term more portable) approach:
 * initialize the global variables.   Only 1 variable in a file needs to
 * be initialized in order for all the variables in the file to be added
 * to the library's index, and thus visible to the linker.
 */
