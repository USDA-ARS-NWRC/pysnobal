#include <errno.h>
#include <math.h>

//#include "ipw.h"
#include "envphys.h"

double
satw(
		double  tk)		/* air temperature (K)		*/
{
	double  x;
	double  l10;

	if (tk <= 0.) {
		fprintf(stderr, "tk=%f\n, less than zero", tk);
//		assert(tk > 0.);
		exit(EXIT_FAILURE);
	}

	errno = 0;
	l10 = log(1.e1);

	x = -7.90298*(BOIL/tk-1.) + 5.02808*log(BOIL/tk)/l10 -
			1.3816e-7*(pow(1.e1,1.1344e1*(1.-tk/BOIL))-1.) +
			8.1328e-3*(pow(1.e1,-3.49149*(BOIL/tk-1.))-1.) +
			log(SEA_LEVEL)/l10;

	x = pow(1.e1,x);

	if (errno) {
		perror("sati: bad return from log or pow");
		exit(EXIT_FAILURE);
//		syserr();
//		error("satw: bad return from log or pow");
	}

	return(x);
}
